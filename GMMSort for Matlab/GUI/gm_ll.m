function [loglike] = gm_ll(x,mu0,S0,alpha0)
    k=length(alpha0);
    n=size(x,1);
    d=size(x,2);
    
    loglike=zeros(k,n);
    for ik=1:k
        try
        L = chol(squeeze(S0(ik,:,:)));

            % L*L' = S0;

            logdetSigma = 2*sum(log(diag(L)));
            % log(det(S0)) = log(det(L)^2) = log(prod(diag(L))^2) = sum(log(diag(L)))*2

    %         x_centered = bsxfun(@minus,x,mu0(ik,:));
            x_centered = x-repmat(mu0(ik,:),[n,1]);


    %         % x-mu
    %         for i=1:n
    %             mahala_2 (i) = x_centered(i,:)*(1./(squeeze(S0(ik,:,:))))*x_centered(i,:)';
    %             % mahalaD = sqrt((x-mu)'*inv(Sigma)*(x-mu))
    %             % mahala_2 = mahalaD^2
    %         end

            xRinv = x_centered/L ;
            mahala_2 = sum(xRinv.^2, 2);

            loglike(ik,:) = -0.5*logdetSigma + log(alpha0(ik)) -0.5*mahala_2 -0.5*d*log(2*pi);
            % log(p(x,z)) = log(alpha0*N(mu0,S0)) = log(alpha0) + log(N(mu0,S0))
        catch e
%             disp(e)
            loglike(ik,:)=nan;
        end
    end
    
end