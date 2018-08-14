function [model, clusterid] = clusterize(features,parameters)
% CLUSTERIZE   Clusterize using the Gaussian mixture models (GMMs)
% and overclustering.
%  
%   [MODEL, CLUSTERID] = CLUSTERIZE(FEATURES,PARAMETERS) where FEATURES is
%   an N-by-M matrix with N samples of M features, and PARAMETERS is a
%   struct containing the fields:
% 
%         ngaussovfit: Number of Gaussians used in the GMM. 
% 
%         nof_replicates: Number of replicates of the model.
% 
%         optgmfit (optional): Struct with the fields max_iter and
%            conv_factor, containing the maximum number of iterations and
%            the threshold for convergence of the GMM, respectively.
% 
%   The function returns:
%
%         MODEL: a struct with the mean (mu), covariance matrix (S) and
%         weight (alpha) of the fitted Gaussians.
% 
%         CLUSTERID: the identity of each sample according to the
%             classification of the model.
% 
% 
% B. C. Souza,
% Brain Institute, Natal, Brazil,
% January, 2018.

%%

if ~isfield(parameters,'optgmfit')
    parameters.optgmfit.max_iter    = 10000;
    parameters.optgmfit.conv_factor = 1e-6;
end


%% Fitting a 5-dim GMM to define the number of clusters and its centers
model=[];
clusterid=[];

nof_dimensions = size(features,2);
flag = 1;
ngaussovfit = parameters.ngaussovfit;
while flag
    try
        MultidimModel = gm_EM(features,ngaussovfit,...
                    'options',parameters.optgmfit,...
                    'replicates',parameters.nof_replicates);
        flag = 0;
    catch
        ngaussovfit = ngaussovfit-1;
        if ngaussovfit <1
            return
        end
    end
end

if size(features,2)==1
    peaks = MultidimModel.mu;
else
    % each Gaussian center is an initial condition to finding peaks
    putativeMultidim_peaks = MultidimModel.mu;
    
    nof_putativepeaks = size(putativeMultidim_peaks,1);
    gradient_minstep = NaN(1,nof_dimensions);
    
    % defining the minimum step for each dimension
    for d = 1 : nof_dimensions
        gradient_minstep(d) = prctile(diff(sort(features(:,d))),1);
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

%% Clustering using the fixed-mean GMM

% Structure telling which parameters are fixed.. in this case, the mean of
% the Gaussians. Can also work fixing S or alpha.
S = struct('mu',peaks,'S',[],'alpha',[]);

model = gm_EM(features,nof_cluster,...
                'options',parameters.optgmfit,...
                'replicates',parameters.nof_replicates,...
                'fixed_param',S);
[~,cluster_probability] = gm_pdf(model,features);

% defining the cluster identity of each spike according to the Gaussian
% with higher probability
[~,clusterid] = max(cluster_probability,[],2);

