function GMM_runSorting( A,B )

global handles
global parameters

%%
waveforms = handles.data.waveforms{handles.chid};
if isfield(handles.dataaux,'nelectrodes')
    nelec = handles.dataaux.nelectrodes;
else
    nelec=1;
end
nof_spikes = size(waveforms,1); 

%% Computing the WD

if(1==1)
    set(handles.Figures.Waveforms.maintext ,'String',...
        'Computing wavelets')
    pause(0.01)
    
    if nelec>1 & mod(size(waveforms,2),nelec)==0
        wavelen = size(waveforms,2)/nelec;
        
        waveforms_WV=[];
        
        for ielec = 1:nelec
            wavetemp = waveforms(:,((ielec-1)*wavelen+1):(ielec*wavelen)); 
        
            N = floor(log2(size(wavetemp,2)));
            [aux, L] = wavedec(wavetemp(1,:),N,'haar');
            wavetemp_WV=nan(size(wavetemp,1), length(aux));
            for i=1:size(wavetemp,1)
                [wavetemp_WV(i,:), L] = wavedec(wavetemp(i,:),N,'haar');
            end
            
            waveforms_WV = cat(2,waveforms_WV,wavetemp_WV);
        end
    else


        N = floor(log2(size(waveforms,2)));    
        [aux, L] = wavedec(waveforms(1,:),N,'haar');
        waveforms_WV=nan(size(waveforms,1), length(aux));
        for i=1:size(waveforms,1)
            [waveforms_WV(i,:), L] = wavedec(waveforms(i,:),N,'haar');
        end
    end
    
    
    % zcoring the wavelet coefficients
    norm_components = zscore(waveforms_WV);
    norm_components(isnan(norm_components)) = 0;

    set(handles.Figures.Waveforms.maintext ,'String',...
        'Computing models for each wavelet: 0%')
    pause(0.01)
    
    nof_wv = size(waveforms_WV,2);
    sep_metric = nan(nof_wv,4);
    param = @(x) [Ipeak(x) Iinf(x)];
    
    gaussMultiFitPDF=nan(nof_wv,1000);
    bins=nan(nof_wv,1000);
    for wavelet_i = 1:nof_wv
        
        componentvalues = norm_components(:,wavelet_i);
        
        ngauss = parameters.maxGauss;
        flag = 1;
        while flag
            try
                
                MultimodalFit = ...
                    gm_EM(componentvalues,ngauss,...
                    'options',parameters.optgmfit,...
                    'replicates',parameters.nof_replicates,...
                    'rep_param',param);
                
                bins(wavelet_i,:) = linspace(min(componentvalues),max(componentvalues),1000);
                gaussMultiFitPDF(wavelet_i,:) = gm_pdf(MultimodalFit,bins(wavelet_i,:)');
                sep_metric(wavelet_i,1:3) = median(MultimodalFit.param);
                
                flag = 0;
            catch errorinfo
                ngauss = ngauss-1;
                if ngauss<1
                    flag = 0;
                    display(errorinfo.message)
                end
            end
        end
        
        set(handles.Figures.Waveforms.maintext ,'String',...
            ['Computing models for each wavelet: ' ...
            num2str(round(wavelet_i*1000/nof_wv)/10) '%'])
        pause(0.01)
        
    end
    
    sep_metric(isnan(sep_metric)) = 0;
    sep_metric(:,4) = var(waveforms_WV);
    
    % % % % % % % % % %
    % choose the separability metric to perform wPCA
    I_peak=1;
    I_inf=2;
    I_dist=3;
    Var=4;
    % % % % % % % % % %
    
    % wPCA is done in the zcored wavelet coefficients
    wcomponents = repmat(sep_metric(:,I_dist)',size(norm_components,1),1).*norm_components;
    wcomponents(isnan(wcomponents)) = 0;
    try
        [~, wpca, ~] = pca(wcomponents);
    catch e
        [~, wpca, ~] = princomp(wcomponents);
    end
    
    %the 5 wPCs are selected for clustering
    Multimodel_components = 1:5;
    waveforms_components = wpca(:,Multimodel_components);
    
    %% PCA and WD
    % % the PCA and WD methods could be done by selecting the 5 best unimodal
    % % components, as folows:
    %[~, idx] = sort(weight,1,'descend');
    %Multimodel_components = idx(1:5,Idist);
    %waveforms_components = waveforms_WV(:,Multimodel_components);
    
    %waveforms = waveforms_WV(:,selected_wavelets);
end


nof_dimensions = length(Multimodel_components);
set(handles.Figures.Waveforms.maintext ,'String',...
    ['Clustering will be perfomed in ' num2str(nof_dimensions) ...
    ' dimension(s)'])
pause(0.01)

%% Fitting a 5-dim GMM to define the number of clusters and its centers

set(handles.Figures.Waveforms.maintext ,'String',...
    'Computing joint PDF'), pause(0.01)
flag = 1;
ngaussovfit = parameters.ngaussovfit;
while flag
    try
        MultidimModel = gm_EM(waveforms_components,ngaussovfit,...
                    'options',parameters.optgmfit,...
                    'replicates',parameters.nof_replicates);
                
        flag = 0;
    catch
        ngaussovfit = ngaussovfit-1;
        if ngaussovfit <1
            set(handles.Figures.Waveforms.maintext ,'String',...
                'Not possible to fit model in multidim space')
            handles.data.clustering_space{handles.chid} = [];
            handles.data.class_id{handles.chid} = nan(1,nof_spikes);
            handles.data.classification_uncertainty{handles.chid} = [];
            return
        end
    end
end

set(handles.Figures.Waveforms.maintext ,'String',...
    'Finding distribution peaks.')
pause(0.01)

if size(waveforms_components,2)==1
    peaks = MultidimModel.mu;
else
    % each Gaussian center is an initial condition to finding peaks
    putativeMultidim_peaks = MultidimModel.mu;
    
    nof_putativepeaks = size(putativeMultidim_peaks,1);
    gradient_minstep = NaN(1,nof_dimensions);
    
    % defining the minimum step for each dimension
    for d = 1 : nof_dimensions
        gradient_minstep(d) = prctile(diff(sort(waveforms_components(:,d))),1);
        % gradient_minstep(d) = 1e-8;
    end
    
    peaks = [];
    
    %%%%%
    % defining the probability of the model as an anonimous function to use
    % fminsearch
    fun_prob_x = @(x) -sum(exp(gm_ll(x,MultidimModel.mu,MultidimModel.S,MultidimModel.alpha)),1);
    %%%%%
    
    for p = 1 : nof_putativepeaks
        % Nelder-Mead algorithm
        current_coord = putativeMultidim_peaks(p,:);
        [x_local_max,fval] = fminsearch(fun_prob_x, current_coord,...
            optimset('TolX',1e-8,'TolFun',1e-8));
        peaks = [peaks; x_local_max];
    end
    
    %normalizing peaks and excluding repeated ones..
    peaks_norm = peaks;
    for d = 1 : nof_dimensions
        peaks_norm(:,d) = round(peaks(:,d)/(5*gradient_minstep(d)));
    end
    [~,centersid]       = unique(peaks_norm,'rows');
    peaks               = peaks(centersid,:);

end

nof_cluster = size(peaks,1);
set(handles.Figures.Waveforms.maintext ,'String',...
    [num2str(nof_cluster) ' clusters detected.'])
pause(0.01)

%% Clustering using the fixed-mean GMM

set(handles.Figures.Waveforms.maintext ,'String',...
    'Clustering.')
pause(0.01)
    
% Structure telling which parameters are fixed.. in this case, the mean of
% the Gaussians. Can also work fixing S or alpha.
S = struct('mu',peaks,'S',[],'alpha',[]);

ClusteringModel = gm_EM(waveforms_components,nof_cluster,...
                'options',parameters.optgmfit,...
                'replicates',parameters.nof_replicates,...
                'fixed_param',S);
[~,cluster_probability] = gm_pdf(ClusteringModel,waveforms_components);

% defining the cluster identity of each spike according to the Gaussian
% with higher probability
[~,id_cluster] = max(cluster_probability,[],2);

if size(cluster_probability,2) ~= nof_cluster
    nof_cluster = size(cluster_probability,2);
end


% computing the entropy of the classes.. this is not used in the interface,
% but can be used to implement a threshold to select only spikes that have
% low overlap of gaussians..
classification_uncertainty = NaN(nof_spikes,1);
for s = 1 : nof_spikes
    classification_uncertainty(s) = 0;
    for c = 1 : nof_cluster
        classification_uncertainty(s) =   classification_uncertainty(s)                                ...
            + (-sum(cluster_probability(s,c)*log2(cluster_probability(s,c) ...
            + ~cluster_probability(s,c))));
    end
end

%%
clusterlabels = unique(id_cluster);

id_cluster_aux = id_cluster;
for clus_i = 1:length(clusterlabels)
    id_cluster_aux(id_cluster==clusterlabels(clus_i)) = clus_i;
end
id_cluster = id_cluster_aux;

% saving the information of the final model, as well as the clustering
% space. This is used for plotting.
handles.data.model{handles.chid} = ClusteringModel;
handles.data.model{handles.chid}.class = 1:length(clusterlabels);
handles.dataaux.model{handles.chid} = handles.data.model{handles.chid};
handles.data.clustering_space{handles.chid} = waveforms_components;
handles.data.class_id{handles.chid} = id_cluster;

handles.data.classification_uncertainty{handles.chid} = ...
    classification_uncertainty;
handles.dataaux.class_id{handles.chid} = id_cluster;


GMM_plotwaveforms
if get(handles.Figures.Waveforms.DispClustersTOGGLE,'value')
    GMM_showclusters
end
figure(handles.Figures.Waveforms.main)

set(handles.Figures.Waveforms.maintext ,'String',...
    ['Sorting done. ' int2str(nof_cluster) ' clusters detected using ' int2str(nof_dimensions) ' componentes.' ])
pause(0.01)

end

