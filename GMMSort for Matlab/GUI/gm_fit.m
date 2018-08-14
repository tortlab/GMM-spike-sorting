function [fit] = gm_fit(gm_model)
% Measure of gaussian separability using the mean bhattacharyya distance
% between the gaussians in the mixture.

k = size(gm_model.mu,1);
if (k>1)
    comb = nchoosek(1:k,2);
    g={};
    for ik=1:k
        g{ik}=struct('mu',gm_model.mu(ik,:),'S',gm_model.S(ik,:,:),'alpha',gm_model.alpha(ik)); 
    end

    fit=zeros(size(comb,1),1);
    for i=1:size(comb,1)
        fit(i)=bhatta(g{comb(i,1)},g{comb(i,2)});
    end
    fit=mean(fit);
else
    fit=0;
end