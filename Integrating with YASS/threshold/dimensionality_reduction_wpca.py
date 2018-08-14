"""
Functions for dimensionality reduction
"""
try:
    from pathlib2 import Path
except Exception:
    from pathlib import Path

from functools import reduce
import logging

import numpy as np
import pywt
from scipy import stats as st
import sklearn.mixture as mix
import os.path 

from yass.batch import BatchProcessor, util
from yass.util import check_for_files, LoadFile, save_numpy_object

logger = logging.getLogger(__name__)


@check_for_files(filenames=[LoadFile('scores_filename'),
                            LoadFile('spike_index_clear_filename'),
                            LoadFile('rotation_matrix_filename')],
                 mode='extract', relative_to='output_path')
def pca(path_to_data, dtype, n_channels, data_order, recordings, spike_index,
        spike_size, temporal_features, neighbors_matrix, channel_index,
        max_memory, gmm_params, output_path=None, scores_filename='scores.npy',
        rotation_matrix_filename='rotation.npy',
        spike_index_clear_filename='spike_index_clear_pca.npy',
        if_file_exists='skip'):
    """Apply PCA in batches

    Parameters
    ----------
    path_to_data: str
        Path to recordings in binary format

    dtype: str
        Recordings dtype

    n_channels: int
        Number of channels in the recordings

    data_order: str
        Recordings order, one of ('channels', 'samples'). In a dataset with k
        observations per channel and j channels: 'channels' means first k
        contiguous observations come from channel 0, then channel 1, and so
        on. 'sample' means first j contiguous data are the first observations
        from all channels, then the second observations from all channels and
        so on

    recordings: np.ndarray (n_observations, n_channels)
        Multi-channel recordings

    spike_index: numpy.ndarray
        A 2D numpy array, first column is spike time, second column is main
        channel (the channel where spike has the biggest amplitude)

    spike_size: int
        Spike size

    temporal_features: numpy.ndarray
        Number of output features

    neighbors_matrix: numpy.ndarray (n_channels, n_channels)
        Boolean numpy 2-D array where a i, j entry is True if i is considered
        neighbor of j

    channel_index: np.array (n_channels, n_neigh)
        Each row indexes its neighboring channels.
        For example, channel_index[c] is the index of
        neighboring channels (including itself)
        If any value is equal to n_channels, it is nothing but
        a space holder in a case that a channel has less than
        n_neigh neighboring channels


    max_memory:
        Max memory to use in each batch (e.g. 100MB, 1GB)

    output_path: str, optional
        Directory to store the scores and rotation matrix, if None, previous
        results on disk are ignored, operations are computed and results
        aren't saved to disk

    scores_filename: str, optional
        File name for rotation matrix if False, does not save data

    rotation_matrix_filename: str, optional
        File name for scores if False, does not save data

    spike_index_clear_filename: str, optional
        File name for spike index clear

    if_file_exists:
        What to do if there is already a file in the rotation matrix and/or
        scores location. One of 'overwrite', 'abort', 'skip'. If 'overwrite'
        it replaces the file if it exists, if 'abort' if raise a ValueError
        exception if the file exists, if 'skip' if skips the operation if the
        file exists

    Returns
    -------
    scores: numpy.ndarray
        Numpy 3D array  of size (n_waveforms, n_reduced_features,
        n_neighboring_channels) Scores for every waveform, second dimension in
        the array is reduced from n_temporal_features to n_reduced_features,
        third dimension depends on the number of  neighboring channels

    rotation_matrix: numpy.ndarray
        3D array (window_size, n_features, n_channels)
    """

    ###########################
    # compute rotation matrix #
    ###########################

    bp = BatchProcessor(path_to_data, dtype, n_channels, data_order,
                        max_memory, buffer_size=spike_size)

    # compute WPCA
    WAVE, FEATURE, CH = 0, 1, 2
    
    logger.info('Preforming WPCA')
    
    logger.info('Computing Wavelets ...')
    feature = bp.multi_channel_apply(wavedec, mode='memory',
                                   spike_index=spike_index,
                                   spike_size=spike_size,
                                   wvtype = 'haar')
    features = reduce(lambda x, y: np.concatenate(x, y), [f for f in feature])
    
    
    logger.info('Computing weights..')
    
    # Weighting the features using metric defined in gmtype
    weights = gmm_weight(features, gmm_params, spike_index)
    wfeatures =  features*weights
    
    n_features = wfeatures.shape[FEATURE]
    wfeatures_lin = np.reshape(wfeatures,(wfeatures.shape[WAVE]*n_features,
                                          wfeatures.shape[CH]))
    feature_index = np.arange(0,wfeatures.shape[WAVE]*n_features,n_features)
    
    
    TMP_FOLDER, _ = os.path.split(path_to_data)
    feature_path = os.path.join(TMP_FOLDER,'features.bin')
    feature_params = writefile(wfeatures_lin, feature_path)
    

    bp_feat = BatchProcessor(feature_path, 
                             feature_params['dtype'], 
                             feature_params['n_channels'], 
                             feature_params['data_order'],
                             max_memory, 
                             buffer_size=n_features)
    
            
    # compute PCA sufficient statistics from extracted features
    
    logger.info('Computing PCA sufficient statistics...')
    stats = bp_feat.multi_channel_apply(suff_stat_features, mode='memory',
                                  spike_index=spike_index,
                                  spike_size=spike_size,
                                  feature_index=feature_index,
                                  feature_size=n_features)
    
         
    suff_stats = reduce(lambda x, y: np.add(x, y), [e[0] for e in stats])
    spikes_per_channel = reduce(lambda x, y: np.add(x, y),
                                [e[1] for e in stats])

    # compute PCA projection matrix
    logger.info('Computing PCA projection matrix...')
    rotation = project(suff_stats, spikes_per_channel, temporal_features,
                       neighbors_matrix)
    
    
    #####################################
    # waveform dimensionality reduction #
    #####################################
    
    
    
    logger.info('Reducing spikes dimensionality with PCA matrix...')
    
    # using a new Batch to read feature file
    res = bp_feat.multi_channel_apply(score_features,
                                 mode='memory',
                                 pass_batch_info=True,
                                 rot=rotation,
                                 channel_index=channel_index,
                                 spike_index=spike_index,
                                 feature_index=feature_index)

    
    scores = np.concatenate([element[0] for element in res], axis=0)
    spike_index = np.concatenate([element[1] for element in res], axis=0)
    feature_index = np.concatenate([element[2] for element in res], axis=0)

    # renormalizing PC projections to similar unitary variance
    scores = st.zscore(scores,axis=0)
    
    # save scores
    if output_path and scores_filename:
        path_to_score = Path(output_path) / scores_filename
        save_numpy_object(scores, path_to_score,
                          if_file_exists=if_file_exists,
                          name='scores')

    if output_path and spike_index_clear_filename:
        path_to_spike_index = Path(output_path) / spike_index_clear_filename
        save_numpy_object(spike_index, path_to_spike_index,
                          if_file_exists=if_file_exists,
                          name='Spike index PCA')

    if output_path and rotation_matrix_filename:
        path_to_rotation = Path(output_path) / rotation_matrix_filename
        save_numpy_object(rotation, path_to_rotation,
                          if_file_exists=if_file_exists,
                          name='rotation matrix')

    return scores, spike_index, rotation

def writefile(var, output_path):
    
    n_channels = var.shape[1]
    f = open(output_path, 'wb')
    var.tofile(f)
    f.close()

    params = util.make_metadata('all', n_channels, str(var.dtype),
                                output_path)
    return params

def wavedec(recordings, spike_index, spike_size, wvtype = 'haar'):
    
    # column ids for index matrix
    SPIKE_TIME, MAIN_CHANNEL = 0, 1
    
    wv = pywt.Wavelet(wvtype)

    window_idx = range(-spike_size, spike_size + 1)
    window_size = len(window_idx)
    
    n_obs, n_channels = recordings.shape
    
    c=0
    
    # computing coefficients for a single waveform
    coeffs = pywt.wavedec(recordings[spike_index[0,SPIKE_TIME]+window_idx,c],
                 wv, mode='constant')
    
    n_features = len(np.hstack(coeffs))
    n_wave = spike_index.shape[0]
    
    # preallocating coeffs
    coeffs = np.zeros((n_wave,n_features,n_channels))

    
    for c in range(n_channels):
        wave = np.zeros((n_wave,window_size))
        for j in range(window_size):
            wave[:,j] = recordings[spike_index[:,SPIKE_TIME]+window_idx[j],c]
        
        coeff = pywt.wavedec(wave, wv, mode='constant')
    
        coeffs[:,:,c] = np.hstack(coeff)
    
    # discarding invalid coefficients
    valid_coeff = np.where(np.var(coeffs[:,:,0],axis=0)!=0);
    
    coeffs=coeffs[:,valid_coeff[0],:]
    
    return coeffs

def suff_stat_features(features, spike_index, spike_size, feature_index, feature_size):
    """
    Get PCA SS matrix of extracted features per recording 
    channel

    Parameters
    ----------
    features: np.ndarray (n_features*n_observations, n_channels)
        Multi-channel extracted features concatenated
    spike_index: np.ndarray (number of spikes, 2)
        Spike indexes as returned from the threshold detector
    spike_size: int
        Spike size
    feature_index: np.array (n_spikes)
        contains the start index of each waveform feature in the 
        linearized feature array
    feature_size: int
        contains the number of extracted features in each waveform 
        (n_features)
        
    Returns
    -------
    numpy.ndarray
        3D array (?)
    numpy.ndarray
        1D array, with n_channels entries, the ith entry contains the number
        of spikes found in the ith channel
    """
    
    
    # column ids for index matrix
    SPIKE_TIME, MAIN_CHANNEL = 0, 1

    n_obs, n_channels = features.shape
    
    window_idx = range(0, feature_size)
    window_size = len(window_idx)
    
    pca_suff_stat = np.zeros((window_size, window_size, n_channels))
    spikes_per_channel = np.zeros(n_channels, 'int32')

    # iterate over every channel
    for c in range(n_channels):
        # get spikes times for the current channel
        channel_spike_times = feature_index[spike_index[:, MAIN_CHANNEL] == c]

        channel_spikes = len(channel_spike_times)

        # create zeros matrix (window size x number of spikes for this channel)
        wf_temp = np.zeros((window_size, channel_spikes))

        # iterate over the window size
        for j in range(window_size):
            # fill in recording values for each spike time
            wf_temp[j, :] = features[channel_spike_times + window_idx[j], c]

        
        pca_suff_stat[:, :, c] = np.matmul(wf_temp, wf_temp.T)

        spikes_per_channel[c] = channel_spikes

    return pca_suff_stat, spikes_per_channel


def suff_stat(recordings, spike_index, spike_size):
    """
    Get PCA SS matrix per recording channel

    Parameters
    ----------
    recordings: np.ndarray (n_observations, n_channels)
        Multi-channel recordings
    spike_index: np.ndarray (number of spikes, 2)
        Spike indexes as returned from the threshold detector
    spike_size: int
        Spike size

    Returns
    -------
    numpy.ndarray
        3D array (?)
    numpy.ndarray
        1D array, with n_channels entries, the ith entry contains the number
        of spikes found in the ith channel
    """
    
    print(recordings.shape)
    # column ids for index matrix
    SPIKE_TIME, MAIN_CHANNEL = 0, 1

    n_obs, n_channels = recordings.shape
    window_idx = range(-spike_size, spike_size + 1)
    window_size = len(window_idx)

    pca_suff_stat = np.zeros((window_size, window_size, n_channels))
    spikes_per_channel = np.zeros(n_channels, 'int32')

    # iterate over every channel
    for c in range(n_channels):
        # get spikes times for the current channel
        channel_spike_times = spike_index[spike_index[:, MAIN_CHANNEL] == c,
                                          SPIKE_TIME]
        channel_spike_times = channel_spike_times[np.logical_and(
            (channel_spike_times > spike_size),
            (channel_spike_times < n_obs - spike_size - 1))]

        channel_spikes = len(channel_spike_times)

        # create zeros matrix (window size x number of spikes for this channel)
        wf_temp = np.zeros((window_size, channel_spikes))

        # iterate over the window size
        for j in range(window_size):
            # fill in recording values for each spike time
            wf_temp[j, :] = recordings[channel_spike_times + window_idx[j], c]

        pca_suff_stat[:, :, c] = np.matmul(wf_temp, wf_temp.T)

        spikes_per_channel[c] = channel_spikes

    return pca_suff_stat, spikes_per_channel

def project(ss, spikes_per_channel, n_features, neighbors):
    """
    Get PCA projection matrix per channel

    Parameters
    ----------
    ss: matrix
        SS matrix as returned from get_pca_suff_stat
    spikes_per_channel: array
        Number of spikes per channel
    n_features: int
        Number of features
    neighbors: matrix
        Neighbors matrix

    Returns
    -------
    numpy.ndarray
        3D array (window_size, n_features, n_channels)
    """
    window_size, _, n_channels = ss.shape
    # allocate rotation matrix for each channel
    rot = np.zeros((window_size, n_features, n_channels))

    ss_all = np.sum(ss, 2)
    w, v = np.linalg.eig(ss_all)
    rot_all = v[:,
                np.argsort(w)[window_size:(window_size - n_features - 1):-1]]

    for c in range(n_channels):
        if spikes_per_channel[c] <= window_size:
            if np.sum(spikes_per_channel[neighbors[c, :]]) <= window_size:
                rot[:, :, c] = rot_all
            else:
                w, v = np.linalg.eig(np.sum(ss[:, :, neighbors[c, :]], 2))
                rot[:, :, c] = v[:,
                                 np.argsort(w)[window_size:(
                                     window_size - n_features - 1):-1]]
        else:
            w, v = np.linalg.eig(ss[:, :, c])
            rot[:, :, c] = v[:,
                             np.argsort(w)[window_size:(
                                 window_size - n_features - 1):-1]]

    return rot


def score(recording, idx_local, idx, rot, channel_index, spike_index):
    """
    Reduce waveform dimensionality using a rotation matrix. Optionally
    return scores only for neighboring channels instead of all channels

    Parameters
    ----------
    recordings: np.ndarray (n_observations, n_channels)
        Multi-channel recordings

    rot: numpy.ndarray
        Rotation matrix. Array with dimensions (n_temporal_features,
        n_features, n_channels) for PCA matrix or (n_temporal_features,
        n_features) for autoencoder matrix

    channel_index: np.array (n_channels, n_neigh)
        Each row indexes its neighboring channels.
        For example, channel_index[c] is the index of
        neighboring channels (including itself)
        If any value is equal to n_channels, it is nothing but
        a space holder in a case that a channel has less than
        n_neigh neighboring channels

    spike_index: np.array (n_spikes, 2)
        contains spike information, the first column is the
        spike time and the second column is the main channel

    Returns
    -------
    scores: np.array (n_spikes, n_features, n_neighboring_channels)
        Scores for every waveform, second dimension in the array is reduced
        from n_temporal_features to n_features, third dimension
        is number of neighboring channels.
    """

    data_start = idx[0].start
    data_end = idx[0].stop
    # get offset that will be applied
    offset = idx_local[0].start

    spike_time = spike_index[:, 0]
    spike_index = spike_index[np.logical_and(spike_time >= data_start,
                                             spike_time < data_end)]
    spike_index[:, 0] = spike_index[:, 0] - data_start + offset

    # obtain shape information
    n_observations, n_channels = recording.shape
    n_data = spike_index.shape[0]
    n_neigh = channel_index.shape[1]

    # if rot has two dimension, rotation matrix is used for every
    # channels, if it is three, the third dimension has to match
    # the number of channels
    if rot.ndim == 2:
        # neural net case
        n_temporal_features, n_reduced_features = rot.shape
        # copy rotation matrix to all channels
        rot = np.tile(rot[:, :, np.newaxis], [1, 1, n_channels])

    elif rot.ndim == 3:
        # pca case
        n_temporal_features, n_features, n_channels_ = rot.shape

        if n_channels != n_channels_:
            raise ValueError('n_channels does not match between '
                             'recording ({}) and the rotation matrix ({})'
                             .format(n_channels,
                                     n_channels_))
    else:
        raise ValueError('rot must have 2 or 3 dimensions (has {})'.format(
            rot.ndim))

    # n_temporal_features has to be an odd number
    if n_temporal_features % 2 != 1:
        raise ValueError('waveform length needs to be'
                         'an odd number (has {})'.format(
                             n_temporal_features))

    R = int((n_temporal_features-1)/2)

    rot = np.transpose(rot, [2, 1, 0])
    scores = np.zeros((n_data, n_features, n_neigh))
    for channel in range(n_channels):

        # get neighboring channel information
        ch_idx = channel_index[channel][
            channel_index[channel] < n_channels]

        # get spikes whose main channel is equal to channel
        idx_c = spike_index[:, 1] == channel

        # get waveforms
        spt_c = spike_index[idx_c, 0]
        waveforms = np.zeros((spt_c.shape[0], ch_idx.shape[0],
                              n_temporal_features))
        for j in range(spt_c.shape[0]):
            waveforms[j] = recording[spt_c[j]-R:spt_c[j]+R+1, ch_idx].T

        # apply rot on wavefomrs
        scores[idx_c, :, :ch_idx.shape[0]] = np.transpose(
            np.matmul(np.expand_dims(rot[ch_idx], 0),
                      np.expand_dims(waveforms, -1))[:, :, :, 0],
            [0, 2, 1])

    spike_index[:, 0] = spike_index[:, 0] + data_start - offset

    return scores, spike_index

def score_features(features, idx_local, idx, rot, channel_index, spike_index, feature_index):
    """
    Similar to score function, but using extracted features instead of 
    standarized recordings. 

    Parameters
    ----------
    features: np.ndarray (n_features*n_spikes, n_channels)
        Multi-channel features of each waveform concatenated.

    rot: numpy.ndarray
        Rotation matrix. Array with dimensions (n_temporal_features,
        n_features, n_channels) for PCA matrix or (n_temporal_features,
        n_features) for autoencoder matrix

    channel_index: np.array (n_channels, n_neigh)
        Each row indexes its neighboring channels.
        For example, channel_index[c] is the index of
        neighboring channels (including itself)
        If any value is equal to n_channels, it is nothing but
        a space holder in a case that a channel has less than
        n_neigh neighboring channels

    spike_index: np.array (n_spikes, 2)
        contains spike information, the first column is the
        spike time and the second column is the main channel
        
    feature_index: np.array (n_spikes)
        contains the start index of each waveform feature in the 
        linearized feature array

    Returns
    -------
    scores: np.array (n_spikes, n_features, n_neighboring_channels)
        Scores for every waveform, second dimension in the array is reduced
        from n_temporal_features to n_features, third dimension
        is number of neighboring channels.
    """

    data_start = idx[0].start
    data_end = idx[0].stop
    # get offset that will be applied
    offset = idx_local[0].start


    valid = np.logical_and(feature_index >= data_start,
                                             feature_index < data_end)
    spike_index = spike_index[valid]
    feature_index =feature_index[valid]
    
    spike_index[:, 0] = spike_index[:, 0] - data_start + offset
    feature_index = feature_index - data_start + offset

    # obtain shape information
    n_observations, n_channels = features.shape
    n_data = feature_index.shape[0]
    n_neigh = channel_index.shape[1]

    # if rot has two dimension, rotation matrix is used for every
    # channels, if it is three, the third dimension has to match
    # the number of channels
    if rot.ndim == 2:
        # neural net case
        n_temporal_features, n_reduced_features = rot.shape
        # copy rotation matrix to all channels
        rot = np.tile(rot[:, :, np.newaxis], [1, 1, n_channels])

    elif rot.ndim == 3:
        # pca case
        n_temporal_features, n_features, n_channels_ = rot.shape

        if n_channels != n_channels_:
            raise ValueError('n_channels does not match between '
                             'recording ({}) and the rotation matrix ({})'
                             .format(n_channels,
                                     n_channels_))
    else:
        raise ValueError('rot must have 2 or 3 dimensions (has {})'.format(
            rot.ndim))

    # n_temporal_features has to be an odd number
    if n_temporal_features % 2 != 1:
        raise ValueError('waveform length needs to be'
                         'an odd number (has {})'.format(
                             n_temporal_features))

    R = int((n_temporal_features-1)/2)

    rot = np.transpose(rot, [2, 1, 0])
    scores = np.zeros((n_data, n_features, n_neigh))
    for channel in range(n_channels):

        # get neighboring channel information
        ch_idx = channel_index[channel][
            channel_index[channel] < n_channels]

        # get spikes whose main channel is equal to channel
        idx_c = spike_index[:, 1] == channel

        # get features
        feat_c = feature_index[idx_c]
        feat = np.zeros((feat_c.shape[0], ch_idx.shape[0],
                              n_temporal_features))
        for j in range(feat_c.shape[0]):
            feat[j] = features[feat_c[j]-R:feat_c[j]+R+1, ch_idx].T

        # apply rot on wavefomrs
        scores[idx_c, :, :ch_idx.shape[0]] = np.transpose(
            np.matmul(np.expand_dims(rot[ch_idx], 0),
                      np.expand_dims(feat, -1))[:, :, :, 0],
            [0, 2, 1])

    spike_index[:, 0] = spike_index[:, 0] + data_start - offset
    feature_index = feature_index + data_start - offset

    
    return scores, spike_index, feature_index

def gmm_weight3(wvcoeff, gmm_params, spike_index):
    
    if 'max_samples' in gmm_params.keys():
        max_n_sample = gmm_params['max_samples']
    else:
        max_n_sample = 10000
        
    if 'replicates' in gmm_params.keys():
        replicates = gmm_params['replicates']
    else:
        replicates = 5
        
    if 'max_iter' in gmm_params.keys():
        max_iter = gmm_params['max_iter']
    else:
        max_iter = 100
        
    if 'n_components' in gmm_params.keys():
        n_components = gmm_params['n_components']
    else:
        n_components = 4
        
    if 'use_channel_features' in gmm_params.keys():
        use_channel_features = gmm_params['use_channel_features']
    else:
        use_channel_features = True
        
    if 'gmtype' in gmm_params.keys():
        gmtype = gmm_params['gmtype']
    else:
        gmtype = 'idist'

    gm = mix.GaussianMixture(covariance_type='full',max_iter=max_iter, n_components=n_components,n_init=1)
    
    
    weight=np.zeros((replicates, wvcoeff.shape[1], wvcoeff.shape[2]))
    for ich in range(wvcoeff.shape[2]):
        
        if use_channel_features:
            idx = np.where(spike_index[:,1]==ich)

            if len(idx[0])>gm.n_components: 
                coeff_channel = wvcoeff[idx]
            else:
                coeff_channel = wvcoeff
        else:
            coeff_channel = wvcoeff
            
        if max_n_sample!=None:
            n_data = coeff_channel.shape[0]
            idx_keep = np.zeros(n_data, 'bool')
            if n_data > max_n_sample:
                idx_sample = np.random.choice(n_data,
                                              max_n_sample,
                                              replace=False)
                idx_keep[idx_sample] = 1
            else:
                idx_keep[:]=1
            coeff_channel=coeff_channel[idx_keep,:,:]

        coeff_norm = st.zscore(coeff_channel,axis=0)
        for irep in range(replicates):
            if gmtype=='idist':
                logger.info('Using distance (idist) metric to weight features..')

                for icoeff in range(coeff_norm.shape[1]):
                    gm.fit(coeff_norm[:,icoeff,ich].reshape(-1,1))
                    dist=[];
                    for i in range(len(gm.means_)):
                        for j in range(i+1,len(gm.means_)):
                            dist.append(np.abs(gm.means_[i]-gm.means_[j])*
                                         np.sqrt(gm.weights_[i]*gm.weights_[j])/
                                         np.sqrt(gm.covariances_[i]*gm.covariances_[j]))

                    weight[irep,icoeff,ich] = np.median(np.hstack(dist))


            elif gmtype=='ipeak':
                logger.info('Using peak (ipeak) metric to weight features..')
                for icoeff in range(coeff_norm.shape[1]):
                    gm.fit(coeff_norm[:,icoeff,ich].reshape(-1,1))
                    model_p = gm.score_samples(np.linspace(np.min(coeff_norm[:,icoeff,ich]),
                                                 np.max(coeff_norm[:,icoeff,ich]),100).reshape(-1,1))

                    idx =detect_peaks(model_p)
                    weight[irep,icoeff,ich] = np.sum(model_p[idx])/np.max(model_p)

            elif gmtype=='iinf':

                logger.info('Using inflection (iinf) metric to weight features..')
                for icoeff in range(coeff_norm.shape[1]):
                    gm.fit(coeff_norm[:,icoeff,ich].reshape(-1,1))
                    model_p = gm.score_samples(np.linspace(np.min(coeff_norm[:,icoeff,ich]),
                                                 np.max(coeff_norm[:,icoeff,ich]),100).reshape(-1,1))
                    #model_p = np.sum(model,axis=1)

                    idx =detect_peaks(np.abs(np.diff(model_p)))
                    weight[irep,icoeff,ich] = np.sum(model_p[idx])/np.max(model_p)

            else:

                logger.info('Incorrect gmtype identifier. Use gmtype= "idist", "ipeak" or "iinf"')

    
    return np.median(weight,axis=0)




def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False, show=False, ax=None):

    """
    Returns
    -------
    ind : 1D array_like
        indeces of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`
    
    The function can handle NaN's 
 
    Copyright 2015. Marcos Duarte. MIT license.

    """

    x = np.atleast_1d(x).astype('float64')
    if x.size < 3:
        return np.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = np.where(np.isnan(x))[0]
    if indnan.size:
        x[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
    ind = np.unique(np.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan-1, indnan+1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size-1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([x[ind]-x[ind-1], x[ind]-x[ind+1]]), axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = np.sort(ind[~idel])
    return ind





















