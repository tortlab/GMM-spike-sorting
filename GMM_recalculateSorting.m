function ClusteringModel = GMM_recalculateSorting (waveforms_components)
global handles
global parameters
%%

nof_dimensions = size(waveforms_components,2);

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
            return
        end
    end
end

%%

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
    end
    
    peaks = [];
    %%%%%
    % defining the probability of the model as an anonimous function to use
    % fminsearch
    fun_prob_x = @(x) -sum(exp(gm_ll(x,MultidimModel.mu,MultidimModel.S,MultidimModel.alpha)),1);
    %%%%%
    tic
    for p = 1 : nof_putativepeaks
        current_coord = putativeMultidim_peaks(p,:);
        [x_local_max,fval] = fminsearch(fun_prob_x, current_coord,...
            optimset('TolX',1e-8,'TolFun',1e-8));
        peaks = [peaks; x_local_max];
    end
    
    peaks_norm = peaks;
    for d = 1 : nof_dimensions
%         peaks_norm(:,d) = round(peaks(:,d)/(5*gradient_minstep));
        peaks_norm(:,d) = round(peaks(:,d)/(5*gradient_minstep(d)));
    end
    [~,centersid]       = unique(peaks_norm,'rows');
    peaks               = peaks(centersid,:);

end

nof_cluster = size(peaks,1);
set(handles.Figures.Waveforms.maintext ,'String',...
    [num2str(nof_cluster) ' clusters detected.'])
pause(0.01)

%% Clustering

set(handles.Figures.Waveforms.maintext ,'String',...
    'Clustering.')
pause(0.01)

ClusteringModel=[];
if nof_cluster >1
    
    S = struct('mu',peaks,'S',[],'alpha',[]);

    
    ClusteringModel = gm_EM(waveforms_components,nof_cluster,...
        'options',parameters.optgmfit,...
        'replicates',parameters.nof_replicates,...
        'fixed_param',S);
    
    if size(ClusteringModel.mu,1) ~= nof_cluster
        nof_cluster = size(ClusteringModel.mu,1);
    end
                
    ClusteringModel.class = 1:nof_cluster;

end

