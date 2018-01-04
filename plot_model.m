function [] = plot_model(model,clusterid,features,dim)
% PLOT_MODEL plots the Gaussians and feature points of the GMM model used
%   to classify. MODEL and CLUSTERID variables are returned by the
%   CLUSTERIZE function, while FEATURES are its input.
%
% 
% B. C. Souza January, 2018
% Brain Institute, Natal, Brazil

dimx=dim(1);
dimy=dim(2);
col = jet(length(model.alpha));
figure();clf

hold on
leg={};
for igauss=1:length(model.alpha)
    plot(features(clusterid==igauss,dimx),features(clusterid==igauss,dimy),'.','color',col(igauss,:))
    leg{igauss} = ['Cluster ' num2str(igauss)];
end

for igauss=1:length(model.alpha)
    
    weigth = 1+10*model.alpha(igauss)...
        /sum(model.alpha);
    offset = model.mu(igauss,[dimx,dimy]);
    aux=linspace(0,2*pi,20);
    x=cos(aux);
    y=sin(aux);
    
    circle = [x;y];
    S = squeeze(model.S(igauss,[dimx,dimy],[dimx,dimy]));
    [U, V] = eig(S);
    elipse = bsxfun(@plus,offset',U*2*sqrt(V)*circle);
    plot(elipse(1,:),elipse(2,:),'color',col(igauss,:),'linewidth',weigth)
    plot(offset(1),offset(2),'o','color','k','markersize',8,'linewidth',4)
    plot(offset(1),offset(2),'o','color',col(igauss,:),'markersize',8,'linewidth',2)
    
end
hold off
axis tight
title('Clustering space', 'fontsize', 16)
xlabel([ 'Feature ' num2str(dimx)], 'fontsize', 14)
ylabel([ 'Feature ' num2str(dimy)], 'fontsize', 14)
legend(leg)
