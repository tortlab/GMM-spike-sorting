function [d_bhatta] = bhatta(g1,g2)
% Bhattacharyya distance between two gaussian distributions

    S = 0.5*g1.S + 0.5*g2.S;
    mu_diff = g1.mu - g2.mu;
    
    L=chol(S);
%     L is actually R=L' that's why we do mu_diff/L instead of mu_diff\L
     mahala2 = sum((mu_diff/L).^2);
    
    logdetSigma = 2*sum(log(diag(L)));
    % log(det(S)) = log(det(L)^2) = log(prod(diag(L))^2) = sum(log(diag(L)))*2
    
    det_S0_S1 = det(g1.S*g2.S);
    % det(s0)*det(s1)= det(s0*s1)
    
    d_bhatta = 0.125*mahala2 + 0.5*logdetSigma - 0.25*log(det_S0_S1);
    

end