function [features] = extract_features(waveforms,parameters,varargin)
% EXTRACT_FEATURES   Extract the features from the waveforms using Gaussian
% mixture models (GMMs) wavelet decomposition and weighted PCA (wPCA).
%  
%   [FEATURES] = EXTRACT_FEATURES(WAVEFORMS,PARAMETERS) where WAVEFORMS is
%   a N-by-M matrix with N waveforms of M points, and PARAMETERS is a
%   struct containing the fields:
% 
%         maxGauss: Number of Gaussians used in the GMM.
% 
%         nof_replicates: Number of replicates of the model.
% 
%         optgmfit (optional): Struct
%            with the fields max_iter and conv_factor, containing the
%            maximum number of iterations and the threshold for convergence
%            of the GMM, respectively.
% 
%   [FEATURES] = EXTRACT_FEATURES(..., OPTION) where OPTION is the
%   separability metric used by the wPCA ('Idist' (default), 'Ipeak',
%   'Iinf' or 'Var')
% 
%   The function returns:
%
%         FEATURES: a matrix containing the wPCA of the wavelet
%         coefficients computed using the separability metric defined.
% 
%  
% B. C. Souza January, 2018
% Brain Institute, Natal, Brazil


% % % % % % % % % %
% choose the separability metric to perform wPCA
used_metric= 3; % Idist;

if ~isempty(varargin)
    if strcmp('Ipeak',varargin{1})
        used_metric=1;
    elseif strcmp('Iinf',varargin{1})
        used_metric=2;
    elseif strcmp('Idist',varargin{1})
        used_metric=3;
    elseif strcmp('Var',varargin{1})
        used_metric=4;
    end
end
% % % % % % % % % %

if ~isfield(parameters,'optgmfit')
    parameters.optgmfit.max_iter    = 10000;
    parameters.optgmfit.conv_factor = 1e-6;
end


N = floor(log2(size(waveforms,2)));
[aux, L] = wavedec(waveforms(1,:),N,'haar');
waveforms_WV=nan(size(waveforms,1), length(aux));
for i=1:size(waveforms,1)
    [waveforms_WV(i,:), L] = wavedec(waveforms(i,:),N,'haar');
end


% zcoring the wavelet coefficients
norm_components = zscore(waveforms_WV);
norm_components(isnan(norm_components)) = 0;

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
            
            UnimodalFit = ...
                gm_EM(componentvalues,ngauss,...
                'options',parameters.optgmfit,...
                'replicates',parameters.nof_replicates,...
                'rep_param',param);
            
            bins(wavelet_i,:) = linspace(min(componentvalues),max(componentvalues),1000);
            gaussMultiFitPDF(wavelet_i,:) = gm_pdf(UnimodalFit,bins(wavelet_i,:)');
            sep_metric(wavelet_i,1:3) = median(UnimodalFit.param,1);
            
            flag = 0;
        catch errorinfo
            ngauss = ngauss-1;
            if ngauss<1
                flag = 0;
%                 display(errorinfo.message)
            end
        end
    end
end

sep_metric(isnan(sep_metric)) = 0;
sep_metric(:,4) = var(waveforms_WV);


% wPCA is done in the zcored wavelet coefficients
wcomponents = repmat(sep_metric(:,used_metric)',size(norm_components,1),1).*norm_components;
wcomponents(isnan(wcomponents)) = 0;
try
    [~, wpca, ~] = pca(wcomponents);
catch e
    [~, wpca, ~] = princomp(wcomponents);
end

features = wpca;
