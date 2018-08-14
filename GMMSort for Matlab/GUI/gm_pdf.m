function [pdf_x, p_z] = gm_pdf(gm_model,x)

k=length(gm_model.alpha);

loglikelihood =gm_ll(x,gm_model.mu,gm_model.S,gm_model.alpha);
pdf_x = sum(exp(loglikelihood),1);


maxll = max(loglikelihood,[],1);
p_z = exp(loglikelihood-repmat(maxll,[k,1]));
% post_prob = p(z_j)*p(x_i|z_j)/exp(maxll)
        
density = sum(p_z,1);
% density = sum_j(p(z_j)*p(x_i|z_j)/exp(maxll)) = p(x_i)/exp(maxll)

p_z = p_z./repmat(density,[k,1]);
p_z = p_z';
% pdf_z = p(z_j)*p(x_i|z_j)/p(x_i) = p(z_j|x_i)

