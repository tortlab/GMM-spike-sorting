function [gm_model] = gm_EM (x, k, varargin)
% Gaussian mixture model using Expectation maximization algorithm
% 
% GM_MODEL = GM_EM(X, K) fits a gaussian mixture model of data in X using K
%   gaussians. X is a N by D matrix, where each row is one observations of
%   a D-dimensional point. GM_MODEL is a struct with the fields:
%           MU: a K by D matrix containing the mean of each gaussian
% 
%           S: a K by D by D matrix containing the K covariance matrices of
%           each gaussian
% 
%           ALPHA: a 1 by K vector containing the weights of each gaussian
% 
%           N_ITER: the number of iteractions done by the algorithm
% 
% GM_MODEL = GM_EM(..., 'options', OPTIONS) specifies OPTIONS, a struct 
%   with the fields
%           MAX_ITER: the max number of iteractions (default = 100)
%           CONV_FACTOR: the convergence factor (default = 1e-6)
% 
% GM_MODEL = GM_EM(..., 'replicates', R) specifies the number of 
%   replications
% 
% GM_MODEL = GM_EM(..., 'fixed_param', PARAM) specifies PARAM, a struct 
%   to fix some parameters of the model. PARAM has fields
%           MU: a K by D matrix as in GM_MODEL. The non NaN values of MU 
%           are fixed in the model.
% 
%           S: a K by D by D matrix as in GM_MODEL. The non NaN values of S 
%           are fixed in the model.
% 
%           ALPHA: a 1 by K vector as in GM_MODEL. The non NaN values of
%           ALPHA are fixed in the model.
% 
% GM_MODEL = GM_EM(..., 'initial_cond', PARAM) specifies PARAM, a struct  
%   to set the initial conditions of the model. PARAM has fields
%           MU: a K by D matrix as in GM_MODEL. 
% 
%           S: a K by D by D matrix as in GM_MODEL.
% 
%           ALPHA: a 1 by K vector as in GM_MODEL. 
% 
%           
% 



%% Initializing the parameters


n=size(x,1);
d=size(x,2);
if isempty(k) 
    k=5;
end

k0 = k;

% set default properties
options.max_iter = 100;
options.conv_factor = 1e-6;

fixed_param.mu =[];
fixed_param.S =[];
fixed_param.alpha =[];


initial_cond.mu =[];
initial_cond.S =[];
initial_cond.alpha =[];

rep =1;

[REG,prop]=parseparams(varargin);


idx=find(strcmpi('Options',prop));
if ~isempty(idx)
    options=prop{idx+1};
end

idx=find(strcmpi('fixed_param',prop));
if ~isempty(idx)
    fixed_param=prop{idx+1};
end

idx=find(strcmpi('initial_cond',prop));
if ~isempty(idx)
    initial_cond=prop{idx+1};
end

idx=find(strcmpi('replicates',prop));
if ~isempty(idx)
    rep=prop{idx+1};
end

param=[];
param_fun=[];
idx=find(strcmpi('rep_param',prop));
if ~isempty(idx)
    param_fun=prop{idx+1};
end

fixed_param0=fixed_param;


%%
gm_model = struct('mu',[],'S',[],'alpha',[],'n_iter',[]);
log_P_x_max=-Inf;
log_P_x_old=-Inf;
log_P_x=-Inf;
for irep = 1:rep
    k=k0;
    fixed_param=fixed_param0;
    
    
% % % % % % % % % % % % % % % % % % %     
%     % using k-means
%     [~, mu0] = kmeans(x,k);
% % % % % % % % % % % % % % % % % % %     
    % set default initial conditions
    mu0 = x(randsample(n,k),:);
    if ~isempty(param_fun)
        mu0=[];
        limits = linspace(min(x),max(x),k+2);
        empty_spaces=0;
        for i=1:length(limits)-2
            idx = find(x>limits(i) & x<= limits(i+2));
            if ~isempty(idx)
                mu0 = [mu0; x(idx(randsample(length(idx),1)))];
            else
                empty_spaces=empty_spaces+1; 
            end
        end
        
        if empty_spaces>0
            mu0 = [mu0; x(randsample(length(x),empty_spaces))];
        end
%         mu0(i) = x(randsample(n,k),:);
    end
% % % % % % % % % % % % % % % % % % % 
    if ~isempty(initial_cond.mu)
        mu0 = initial_cond.mu;
    end
    
    S0=zeros(k,d,d);
    for ik=1:k
        S0(ik,:,:) = cov(x);
    end
    if ~isempty(initial_cond.S)
        S0 = initial_cond.S;
    end
    
    alpha0=1/k*ones(k,1);
    if ~isempty(initial_cond.alpha)
        alpha0 = initial_cond.alpha;
    end

    for iiter = 1:options.max_iter
    %% Computing the Log likelihood function

        loglikelihood = gm_ll(x,mu0,S0,alpha0);
        
        % Exclude any component ik with probability lower than 0.001% 
        if sum(alpha0<0.00001) || sum(isnan(loglikelihood(:,1)))
            if sum(isnan(loglikelihood(:,1)))
                idx = find(~isnan(loglikelihood(:,1)));
            else
                idx=find(alpha0>=0.00001);
            end
            alpha0 = alpha0(idx);
            mu0=mu0(idx,:);
            S0=S0(idx,:,:);
            loglikelihood=loglikelihood(idx,:);
            
            if ~isempty(fixed_param.alpha)
                fixed_param.alpha=fixed_param.alpha(idx);
            end
            if ~isempty(fixed_param.mu)
                fixed_param.mu=fixed_param.mu(idx,:);
            end
            if ~isempty(fixed_param.S)
                fixed_param.S=fixed_param.S(idx,:,:);
            end
            
            % update number of components
            k=length(idx);
        end
        
        if k==0
            break;
        end
        %% E-step
        
        maxll = max(loglikelihood,[],1);

        post_prob = exp(loglikelihood-repmat(maxll,[k,1]));
%         post_prob = exp(bsxfun(@minus,loglikelihood,maxll));        
        % post_prob = p(z_j)*p(x_i|z_j)/exp(maxll)
        
        density = sum(post_prob,1);
        % density = sum_j(p(z_j)*p(x_i|z_j)/exp(maxll)) = p(x_i)/exp(maxll)

        post_prob = post_prob./repmat(density,[k,1]);
%         post_prob = bsxfun(@rdivide,post_prob,density);
        % post_prob = p(z_j)*p(x_i|z_j)/p(x_i) = p(z_j|x_i)

        log_P_x = sum(log(density) +maxll);
        % log_prob_x = log(p(x_i)/exp(maxll)) + maxll = log(p(x_i))-maxll+maxll=
        % = log(p(x_i))

        
        %% Check convergence
        
        log_P_diff = log_P_x-log_P_x_old;
        if log_P_diff>=0 & log_P_diff < options.conv_factor *abs(log_P_x)
            break;
        end
        log_P_x_old = log_P_x;
        %% M-step

        n_k = sum(post_prob,2);
        alpha = n_k/sum(n_k);
        

        mu=zeros(k,d);
        S=zeros(k,d,d);

        for ik=1:k
            mu(ik,:) = post_prob(ik,:)*x/n_k(ik);

            x_centered = x - repmat(mu(ik,:),[n,1]);
            x_centered = repmat(sqrt(post_prob(ik,:)),[d,1])'.*x_centered;
%             x_centered = bsxfun(@minus,x,mu(ik,:));
%             x_centered = bsxfun(@times,sqrt(post_prob(ik,:))', x_centered);

            S(ik,:,:) = x_centered'*x_centered/n_k(ik) *eye(d);

        %     temp =zeros(d,d);
        %     for i=1:n
        %         temp = temp + post_prob(ik,i)*x(i,:)'*x(i,:);
        %     end
        %     S2(ik,:,:) = temp/n_k(ik)-mu(ik,:)'*mu(ik,:) *eye(d);

        end
        %% Actualizing the values

        alpha0=alpha;
        mu0=mu;
        S0=S;

        %% Checking if there are fixed parameters

        if ~isempty(fixed_param.mu)
            mu0(~isnan(fixed_param.mu)) = fixed_param.mu(~isnan(fixed_param.mu));
        end

        if ~isempty(fixed_param.S)
            S0(~isnan(fixed_param.S)) = fixed_param.S(~isnan(fixed_param.S));
        end

        if ~isempty(fixed_param.alpha)
            alpha0(~isnan(fixed_param.alpha)) = fixed_param.alpha(~isnan(fixed_param.alpha));

            % Corrects alpha0 so that sum(alpha0)=1
            scale_factor = (1-sum(alpha0(~isnan(fixed_param.alpha))))/sum(alpha0(isnan(fixed_param.alpha)));
            alpha0(isnan(fixed_param.alpha)) = alpha0(isnan(fixed_param.alpha))*scale_factor;

            sum(alpha0)
        end
        
    end
    %%
    
    if log_P_x>log_P_x_max
        gm_model.mu = mu0;
        gm_model.S = S0;
        gm_model.alpha = alpha0;
        gm_model.n_iter = iiter;
        log_P_x_max=log_P_x;
    end
    
    temp_model.mu = mu0;
    temp_model.S = S0;
    temp_model.alpha = alpha0;
    temp_model.n_iter = iiter;
    
    distance=[];
    len = min([length(temp_model.mu) length(temp_model.alpha)]);
    for igauss = 1:(len-1)
        for jgauss = (igauss+1):len
            distance = [distance sqrt(temp_model.alpha(igauss)*temp_model.alpha(jgauss))*...
                abs(temp_model.mu(igauss)-temp_model.mu(jgauss))...
                /sqrt(temp_model.S(igauss)*temp_model.S(jgauss))];
        end
    end    
    Idist=median(distance);
    
    if ~isempty(param_fun)
        bins=linspace(min(x),max(x),1000);
        pdf = gm_pdf(temp_model,bins');          
        param(irep,:) = [param_fun(pdf) Idist ];
    end
    
end

gm_model.param = param;







%%