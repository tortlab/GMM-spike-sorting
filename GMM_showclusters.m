function GMM_showclusters( A,B )

global handles

if 1==1 %isfield(handles.data,'clustering_space')
    if ~isfield(handles.Figures,'Clusters') || ~isfield(handles.Figures.Clusters,'main') 
        handles.Figures.Clusters.main = figure('units','normalized','position',[.1 .1 .7 .8]);clf
        set(handles.Figures.Clusters.main, 'MenuBar', 'none');

        GMM_drawclusters
    else
        GMM_resetclusters
    end
    GMM_plotclusteringspace(1)

% %  Use this to always activate the stability and crosscorrelation
% %  analyses 
%     GMM_plotstability(3)
%     GMM_runxcorr

end
end

function GMM_drawclusters( A,B )
global handles

figure(handles.Figures.Clusters.main)
set(gcf,'color','w')
positionaux = [0.04 .45 .4 .45];
handles.Figures.Clusters.subplots(1) = ...
    subplot('Position',positionaux);cla

positionaux = [0.55 .45 .4 .45];
handles.Figures.Clusters.subplots(2) = ...
    subplot('Position',positionaux);cla


handles.Figures.Clusters.currentaxes = [1 2];


nof_dimensions = size(handles.data.clustering_space{handles.chid},2);
if nof_dimensions==1
    handles.Figures.Clusters.currentaxes = [1 1];
end

subplotaux = 1;
positionaux = [0.005 .525 .023 .375];
handles.Figures.Clusters.display_y(subplotaux) = ...
    uicontrol('Style', 'popup',...
    'Units','normalized',...
    'String', 1:nof_dimensions,...
    'Position', positionaux,...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'Callback',@GMM_set_xydim_cluster,...
    'fontsize',15,...
    'horizontalAlignment', 'center');
set(handles.Figures.Clusters.display_y(subplotaux),'Value',...
    handles.Figures.Clusters.currentaxes(1))

positionaux = [0.42 .045 .023 .375];
handles.Figures.Clusters.display_x(subplotaux) = ...
    uicontrol('Style', 'popup',...
    'Units','normalized',...
    'String', 1:nof_dimensions,...
    'Position', positionaux,...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'Callback',@GMM_set_xydim_cluster,...
    'fontsize',15,...
    'horizontalAlignment', 'center');
set(handles.Figures.Clusters.display_x(subplotaux),'Value',...
    handles.Figures.Clusters.currentaxes(2))



classlabels = handles.data.class_id{handles.chid};
classlabels = unique(classlabels(~isnan(classlabels)));

popupauxGSF = cellfun(@num2str,num2cell(classlabels),'Uniformoutput',0);
popupauxGSF{end+1} = 'all';
nclasses = length(classlabels);

handles.Figures.Clusters.getspikesfrom = uicontrol('units','normalized',...
    'Style','popup',...
    'Position',[0.17 0.93 0.045 0.02],...
    'String',popupauxGSF,...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'fontsize',15);
set(handles.Figures.Clusters.getspikesfrom,'value',nclasses+1)


handles.Figures.Clusters.viewspikesText = uicontrol('units','normalized',...
    'Style','text',...
    'Position',[0.45 0.83 0.045 0.04],...
    'String','Plot:',...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'fontsize',15);

handles.Figures.Clusters.viewspikesList = uicontrol('units','normalized',...
    'Style','list','max',8,'min',1,...
    'Position',[0.45 0.69 0.05 0.14],...
    'String',popupauxGSF,...
    'Value',1,...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'Callback',@GMM_ChangeView,...
    'fontsize',15);


if isfield(handles.Figures.Clusters,'mode_status')
    switch handles.Figures.Clusters.mode_status
        case 1
            mode_aux='Waveforms';
        case 2
            mode_aux='Components';
        case 3
            mode_aux='Clusters';
    end
else
    if nclasses>0
        set(handles.Figures.Clusters.viewspikesList,'value',nclasses+1)
    end
    handles.Figures.Clusters.mode_status = 1;
    mode_aux = 'Waveforms';
end

handles.Figures.Clusters.changemodeBUTTON = uicontrol('units','normalized',...
    'Style','pushbutton',...
    'Position',[0.27 0.92 0.1 0.04],...
    'String',mode_aux,...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'Callback',@GMM_ChangeMode,...
    'fontsize',14,...
    'horizontalAlignment', 'center');

handles.Figures.Clusters.changevisualizationBUTTON = uicontrol('units','normalized',...
    'Style','togglebutton',...
    'Position',[0.38 0.92 0.05 0.04],...
    'String','%',...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'Callback',@GMM_ChangeVisualization,...
    'fontsize',14,...
    'horizontalAlignment', 'center');



if isfield(handles.Figures.Clusters,'viewspikesVector')
    set(handles.Figures.Clusters.viewspikesList,'value',handles.Figures.Clusters.viewspikesVector)
else
    set(handles.Figures.Clusters.viewspikesList,'value',nclasses+1)
    handles.Figures.Clusters.viewspikesVector = classlabels;
end


handles.Figures.Clusters.mergeButton = uicontrol('units','normalized',...
    'Style','pushbutton',...
    'Position',[0.45 0.64 0.05 0.04],...
    'String','Merge',...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'Callback',@GMM_mergeclusters,...
    'fontsize',15);

handles.Figures.Clusters.deleteButton = uicontrol('units','normalized',...
    'Style','pushbutton',...
    'Position',[0.45 0.58 0.05 0.04],...
    'String','Delete',...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'Callback',@GMM_deleteclusters,...
    'fontsize',15);

handles.Figures.Clusters.recalculateButton = uicontrol('units','normalized',...
    'Style','pushbutton',...
    'Position',[0.45 0.52 0.05 0.04],...
    'String','ReSort',...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'Callback',@GMM_recalculateclusters,...
    'fontsize',15);


handles.Figures.Clusters.moveto = uicontrol('units','normalized',...
    'Style','pushbutton',...
    'Position',[0.55 0.92 0.09 0.04],...
    'Callback', @GMM_MoveTo,...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'String','Move spikes to:',...
    'fontsize',15);

classlabels = handles.data.class_id{handles.chid};
classlabels = unique(classlabels(~isnan(classlabels)));

popupauxBCK = {num2str(classlabels)};
popupauxBCK{end+1} = 'unsorted';
popupauxBCK{end+1} = 'new cluster';
% popupauxBCK{end+1} = '';
nclasses = length(classlabels);

handles.Figures.Clusters.selectedBKGNtemplateH = uicontrol('units','normalized',...
    'Style','popup',...
    'Position',[0.65 0.93 0.09 0.02],...
    'String',popupauxBCK,...
    'callback',@GMM_change_bck_template,...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'fontsize',15); 
set(handles.Figures.Clusters.selectedBKGNtemplateH,'value',nclasses+2)

GMM_setselectedspks(zeros(length(handles.data.class_id{handles.chid}),1)==1);

handles.Figures.Clusters.startpolygon = ...
    uicontrol ('Style', 'pushbutton', 'String', 'Select from:','Units',...
    'normalized','Position',[0.075 0.92 0.085 0.04],...
    'Callback', @GMM_startpolygon,...
    'fontsize',14,...
    'BackgroundColor', [1 1 1]*.9,'ForegroundColor', [1 1 1]*0);


positionaux = [0.04 .05 .4 .3];
handles.Figures.Clusters.crosscorrelation = ...
    subplot('Position',positionaux);
axes(handles.Figures.Clusters.crosscorrelation)
ylabel('r','fontsize',16)
xlabel('Time (ms)','fontsize',16)
title('Cross-correlation','fontsize',18)

handles.Figures.Clusters.crosscorrelationN1 = uicontrol('units','normalized',...
    'Style','popup',...
    'Position',[0.005 0.335 0.023 0.02],...
    'String',{num2str(classlabels)},...
    'Value',1,...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'fontsize',15,...
    'callback',@GMM_runxcorr);

handles.Figures.Clusters.crosscorrelationN2 = uicontrol('units','normalized',...
    'Style','popup',...
    'Position',[0.005 0.3 0.023 0.02],...
    'String',{num2str(classlabels)},...
    'Value',1,...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'fontsize',15,...
    'callback',@GMM_runxcorr);



positionaux = [0.55 .05 .4 .3];
handles.Figures.Clusters.subplots(3) = ...
    subplot('Position',positionaux);cla
ylabel('Peak-to-valley','fontsize',16)
xlabel('Time (min)','fontsize',16)
title('Cluster Stability','fontsize',18)

ax2 = axes('Position', get(handles.Figures.Clusters.subplots(3),'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Xtick',[],...
    'Color','none');

handles.Figures.Clusters.subplots(4) = ax2;


handles.Figures.Clusters.stability = uicontrol('units','normalized',...
    'Style','popup',...
    'Position',[0.5065 0.335 0.023 0.02],...
    'String',{num2str(classlabels)},...
    'Value',1,...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'fontsize',15,...
    'callback',@GMM_changestability);

end


function GMM_updateclusters( A,B )
global handles

figure(handles.Figures.Clusters.main)

axes(handles.Figures.Clusters.subplots(1));cla
axes(handles.Figures.Clusters.subplots(2));cla

if ~isfield(handles.Figures.Clusters,'currentaxes')
    handles.Figures.Clusters.currentaxes = [1 2];
end

nof_dimensions = size(handles.data.clustering_space{handles.chid},2);

subplotaux = 1;


if sum(handles.Figures.Clusters.currentaxes(1:2)>nof_dimensions)
    handles.Figures.Clusters.currentaxes(1:2)=1;
end
set(handles.Figures.Clusters.display_y(subplotaux),...
    'String', 1:nof_dimensions,...
    'Value', handles.Figures.Clusters.currentaxes(1))

set(handles.Figures.Clusters.display_x(subplotaux),...
    'String', 1:nof_dimensions,...
    'Value', handles.Figures.Clusters.currentaxes(2))


classlabels = handles.data.class_id{handles.chid};
classlabels = unique(classlabels(~isnan(classlabels)));
popupauxGSF = {num2str(classlabels)};
popupauxGSF{end+1} = 'all';
nclasses = length(classlabels);

set(handles.Figures.Clusters.getspikesfrom,...
    'String',popupauxGSF,...
    'value',nclasses+1)


set(handles.Figures.Clusters.viewspikesList,...
    'String',popupauxGSF)

set(handles.Figures.Clusters.viewspikesList,...
    'Value',length(handles.Figures.Clusters.viewspikesList))


if isfield(handles.Figures.Clusters,'mode_status')
    switch handles.Figures.Clusters.mode_status
        case 1
            mode_aux='Waveforms';
        case 2
            mode_aux='Components';
        case 3
            mode_aux='Clusters';
    end
else
    set(handles.Figures.Clusters.viewspikesList,'value',nclasses+1)
    handles.Figures.Clusters.mode_status = 1;
    mode_aux = 'Waveforms';
end

set(handles.Figures.Clusters.changemodeBUTTON,...
    'String',mode_aux)
    
if isfield(handles.Figures.Clusters,'viewspikesVector')
    if(~isnan(handles.Figures.Clusters.viewspikesVector))
        set(handles.Figures.Clusters.viewspikesList,'value',handles.Figures.Clusters.viewspikesVector)
    end
else
    set(handles.Figures.Clusters.viewspikesList,'value',nclasses+1)
    handles.Figures.Clusters.viewspikesVector = classlabels;
end


popupauxBCK = {num2str(classlabels)};
popupauxBCK{end+1} = 'unsorted';
popupauxBCK{end+1} = 'new cluster';

set(handles.Figures.Clusters.selectedBKGNtemplateH,...
    'String',popupauxBCK,...
    'value',nclasses+2)

GMM_setselectedspks(zeros(length(handles.data.class_id{handles.chid}),1)==1);

axes(handles.Figures.Clusters.crosscorrelation)
cla

axes(handles.Figures.Clusters.subplots(3))
cla
axes(handles.Figures.Clusters.subplots(4))
cla

set(handles.Figures.Clusters.crosscorrelationN1,...
    'String',{num2str(classlabels)},...
    'Value',1);

set(handles.Figures.Clusters.crosscorrelationN2,...
    'String',{num2str(classlabels)},...
    'Value',1);


set(handles.Figures.Clusters.stability,...
    'String',{num2str(classlabels)},...
    'Value',1);

end


function GMM_resetclusters( A,B )
global handles


classlabels = handles.data.class_id{handles.chid};
classlabels = unique(classlabels(~isnan(classlabels)));
nclasses = length(classlabels);

% set(handles.Figures.Clusters.viewspikesList,'value',nclasses+1)
handles.Figures.Clusters.viewspikesVector = classlabels;

GMM_updateclusters

end

function [] = GMM_recalculateclusters(A,B)
    
global handles

    
clusters = handles.Figures.Clusters.viewspikesVector;
clusters = sort(clusters);

clusterid = handles.data.class_id{handles.chid};
idx = 0*clusterid==1;

idel= handles.data.model{handles.chid}.class*0==1;

for i_class=1:length(clusters)
    idx=idx | clusterid==clusters(i_class);
    
    idel= idel | handles.data.model{handles.chid}.class==clusters(i_class);

end

new_model = GMM_recalculateSorting(handles.data.clustering_space{handles.chid}(idx,:));

if ~isempty (new_model)

    handles.data.model{handles.chid}.mu=handles.data.model{handles.chid}.mu(~idel,:);
    handles.data.model{handles.chid}.alpha=handles.data.model{handles.chid}.alpha(~idel);
    handles.data.model{handles.chid}.S=handles.data.model{handles.chid}.S(~idel,:,:);
    handles.data.model{handles.chid}.class=handles.data.model{handles.chid}.class(~idel);
    nclass = max(handles.data.model{handles.chid}.class);
    if isempty(nclass)
        nclass=0;
    end
    
    handles.data.model{handles.chid}.mu = cat(1,handles.data.model{handles.chid}.mu,new_model.mu);
    handles.data.model{handles.chid}.alpha = cat(1,handles.data.model{handles.chid}.alpha,...
        new_model.alpha*(1-sum(handles.data.model{handles.chid}.alpha)));
    handles.data.model{handles.chid}.S = cat(1,handles.data.model{handles.chid}.S,new_model.S);
    handles.data.model{handles.chid}.class = [handles.data.model{handles.chid}.class...
        nclass+new_model.class];

    [~,cluster_probability] = gm_pdf(handles.data.model{handles.chid},handles.data.clustering_space{handles.chid}(:,:));
    [~,id_cluster] = max(cluster_probability,[],2);

    clusterid(:) = handles.data.model{handles.chid}.class(id_cluster);
    
    clusterlabels = unique(clusterid(~isnan(clusterid)));

    for i = 1:length(clusterlabels);

        clusspikes = clusterid==clusterlabels(i);
        handles.data.class_id{handles.chid}(clusspikes) = i;

        handles.data.model{handles.chid}.class(handles.data.model{handles.chid}.class==clusterlabels(i)) = i;

    end
    handles.Figures.Clusters.viewspikesVector = 1:length(clusterlabels);

    %%

    GMM_plotwaveforms
    if get(handles.Figures.Waveforms.DispClustersTOGGLE,'value')
        GMM_showclusters
    end

end
end

function [ ] = GMM_deleteclusters(A,B)

global handles


    %%
    
clusters = handles.Figures.Clusters.viewspikesVector;
clusters = sort(clusters);

clusterid = handles.data.class_id{handles.chid};
idx = 0*clusterid==1;

idel= handles.data.model{handles.chid}.class*0==1;

for i_class=1:length(clusters)
    idx=idx | clusterid==clusters(i_class);
    
    idel= idel | handles.data.model{handles.chid}.class==clusters(i_class);

end

handles.data.model{handles.chid}.mu=handles.data.model{handles.chid}.mu(~idel,:);
handles.data.model{handles.chid}.alpha=handles.data.model{handles.chid}.alpha(~idel)/...
    sum(handles.data.model{handles.chid}.alpha(~idel));
handles.data.model{handles.chid}.S=handles.data.model{handles.chid}.S(~idel,:,:);
handles.data.model{handles.chid}.class=handles.data.model{handles.chid}.class(~idel);



[~,cluster_probability] = gm_pdf(handles.data.model{handles.chid},handles.data.clustering_space{handles.chid}(idx,:));
% [~,cluster_probability(:,~idel)] = gm_pdf(handles.data.model{handles.chid},handles.data.clustering_space{handles.chid}(idx,:));
[~,id_cluster] = max(cluster_probability,[],2);

clusterid(idx) = handles.data.model{handles.chid}.class(id_cluster);

% handles.data.model{handles.chid}.class(clusters) = clusters(1);

clusterlabels = unique(clusterid(~isnan(clusterid)));

for i = 1:length(clusterlabels);
    
    clusspikes = clusterid==clusterlabels(i);
    handles.data.class_id{handles.chid}(clusspikes) = i;
    
    handles.data.model{handles.chid}.class(handles.data.model{handles.chid}.class==clusterlabels(i)) = i;
    
end
handles.Figures.Clusters.viewspikesVector = 1:length(clusterlabels);

%%


GMM_plotwaveforms
if get(handles.Figures.Waveforms.DispClustersTOGGLE,'value')
    GMM_showclusters
end

end


function [ ] = GMM_mergeclusters(A,B)

global handles

clusters = handles.Figures.Clusters.viewspikesVector;
clusters = sort(clusters);

clusterid = handles.data.class_id{handles.chid};
idx = 0*clusterid==1;
for i=1:length(clusters)
    final_model = min(find(handles.data.model{handles.chid}.class==clusters(i)));
    if ~isempty(final_model)
        break;
    end
end
for i_class=1:length(clusters)
    idx=idx | clusterid==clusters(i_class);
    
    
    handles.data.model{handles.chid}.class(handles.data.model{handles.chid}.class==clusters(i_class)) = ...
        handles.data.model{handles.chid}.class(final_model);
end
clusterid(idx)=clusters(i);


clusterlabels = unique(clusterid(~isnan(clusterid)));

for i = 1:length(clusterlabels);
    
    clusspikes = clusterid==clusterlabels(i);
    handles.data.class_id{handles.chid}(clusspikes) = i;
    
    handles.data.model{handles.chid}.class(handles.data.model{handles.chid}.class==clusterlabels(i)) = i;
    
end
handles.Figures.Clusters.viewspikesVector = 1:length(clusterlabels);

GMM_plotwaveforms
if get(handles.Figures.Waveforms.DispClustersTOGGLE,'value')
    GMM_showclusters
end


end

function [] = GMM_CheckActiveElements()
global handles

set(handles.Figures.Clusters.display_x(1),'enable','on');
set(handles.Figures.Clusters.display_y(1),'enable','on');
set(handles.Figures.Clusters.getspikesfrom,'enable', 'on');
set(handles.Figures.Clusters.changevisualizationBUTTON,'enable', 'on');
set(handles.Figures.Clusters.startpolygon,'enable','on')

if (get(handles.Figures.Clusters.changevisualizationBUTTON,'value'))
    set(handles.Figures.Clusters.startpolygon,'enable','off')
end

switch handles.Figures.Clusters.mode_status
    case 1
        % Showing Clusters!
        set(handles.Figures.Clusters.changevisualizationBUTTON,'enable', 'on');
        set(handles.Figures.Clusters.startpolygon,'enable','on')
    case 2
        % Showing Waveforms!
        set(handles.Figures.Clusters.getspikesfrom,'enable', 'off');
        set(handles.Figures.Clusters.display_x(1),'enable','off');
        set(handles.Figures.Clusters.display_y(1),'enable','off');
        
        
    case 3
        % Showing Components!
        set(handles.Figures.Clusters.getspikesfrom,'enable', 'off');
        set(handles.Figures.Clusters.display_x(1),'enable','off');
        set(handles.Figures.Clusters.display_y(1),'enable','off');
        
        
end



end

function [] = GMM_ChangeVisualization(A,B)
    global handles

    if get(handles.Figures.Clusters.changevisualizationBUTTON,'value')
        set(handles.Figures.Clusters.changevisualizationBUTTON,'string','=')
    else
        set(handles.Figures.Clusters.changevisualizationBUTTON,'string','%')
    end
    

    GMM_plotclusteringspace(1);

end

function [] = GMM_ChangeMode(A,B)

global handles


switch handles.Figures.Clusters.mode_status
    case 1
        % Show Waveforms!
        handles.Figures.Clusters.mode_status=2;
        next_mode='Components';
        
    case 2
        % Show Components!
        handles.Figures.Clusters.mode_status=3;
        next_mode='Clusters';

    case 3
       
        % Show Clusters!
        handles.Figures.Clusters.mode_status=1;
        next_mode='Waveforms';
        
end

set(A,'String',next_mode)


GMM_plotclusteringspace(1)

end

function [] = GMM_ChangeView(A,B)

global handles

classlabels = handles.data.class_id{handles.chid};
classlabels = unique(classlabels(~isnan(classlabels)));
nclasses = length(classlabels);
    
if get(A,'value')>nclasses
handles.Figures.Clusters.viewspikesVector=1:nclasses;
else
handles.Figures.Clusters.viewspikesVector=get(A,'value');
end

GMM_plotclusteringspace(1)

end

function [] = GMM_set_xydim_cluster(A,B)

global handles

subplotaux = [find(handles.Figures.Clusters.display_y==A) ...
    find(handles.Figures.Clusters.display_x==A)];

GMM_plotclusteringspace(subplotaux)


end

function [] = GMM_plotclusteringspace(subplotaux)

global handles


GMM_CheckActiveElements

axes(handles.Figures.Clusters.subplots(subplotaux))
cla 


dimx = get(handles.Figures.Clusters.display_x(subplotaux),'value');
dimy = get(handles.Figures.Clusters.display_y(subplotaux),'value');


handles.Figures.Clusters.currentaxes = [dimy dimx];

class_id = handles.data.class_id{handles.chid};
class_id = class_id(~isnan(class_id));
classlabels = unique(class_id);
nclasses = length(classlabels);
if sum(~isnan(handles.data.clustering_space{handles.chid}(1:10)))==0
    handles.Figures.Clusters.mode_status =2;
    set(handles.Figures.Clusters.changemodeBUTTON,'string','Components');
end
hold on
for class_i = 1:nclasses
%     handles.Figures.Clusters.viewspikesVector
%     classlabels
    if ismember(class_i,handles.Figures.Clusters.viewspikesVector)
        plotspikes = handles.data.class_id{handles.chid}==class_i;
        auxplot = find(plotspikes);
        wave_max = 1000;
        nwave = length(auxplot);
        if handles.Config.Plot==1
            
            if nwave>wave_max
                [v, irand]=sort(rand(1,length(auxplot)));
                auxplot = auxplot(irand(1:wave_max));

            end            
        end
        
        switch handles.Figures.Clusters.mode_status
            case 1
                if(max([dimx dimy])>size(handles.data.clustering_space{handles.chid},2))
                    set(handles.Figures.Clusters.display_x(subplotaux),'value',1);
                    set(handles.Figures.Clusters.display_y(subplotaux),'value',1);
                    dimx=1;
                    dimy=1;
                end

                if dimx~=dimy
                    plot(handles.data.clustering_space{handles.chid}(auxplot,dimx),...
                        handles.data.clustering_space{handles.chid}(auxplot,dimy),'.',...
                        'color',handles.dataaux.class_colors(class_i,:))
                    
                    try
                        for i_model=find(handles.data.model{handles.chid}.class==class_i)


                            weigth = 1+10*handles.data.model{handles.chid}.alpha(i_model)...
                                /sum(handles.data.model{handles.chid}.alpha);
                            offset = handles.data.model{handles.chid}.mu(i_model,[dimx,dimy]);
                            aux=linspace(0,2*pi,20);
                            x=cos(aux);
                            y=sin(aux);

                            circle = [x;y];
                            S = squeeze(handles.data.model{handles.chid}.S(i_model,[dimx,dimy],[dimx,dimy]));
                            [U, V] = eig(S);
                            elipse = bsxfun(@plus,offset',U*2*sqrt(V)*circle);
                            plot(elipse(1,:),elipse(2,:),'color',handles.dataaux.class_colors(class_i,:),'linewidth',weigth)
                            plot(offset(1),offset(2),'o','color','k','markersize',8,'linewidth',4)
                            plot(offset(1),offset(2),'o','color',handles.dataaux.class_colors(class_i,:),'markersize',8,'linewidth',2)
                        end
                    catch e
                        disp('Impossible to plot model elipses.')
                    end
                    
%                     xlim(minmax(handles.data.clustering_space{handles.chid}(:,dimx)'))
%                     ylim(minmax(handles.data.clustering_space{handles.chid}(:,dimy)'))
                else
                    [counts edges] = hist(handles.data.clustering_space{handles.chid}(auxplot,dimx),50);
                    stem(edges(counts~=0),counts(counts~=0),'color',handles.dataaux.class_colors(class_i,:))
                end
                
            case 2
                if get(handles.Figures.Clusters.changevisualizationBUTTON,'value')
                    GMM_PlotPrctiles(handles.data.waveforms{handles.chid}(plotspikes,:)',...
                        handles.dataaux.class_colors(class_i,:))
                else
                    plot(handles.data.waveforms{handles.chid}(auxplot,:)',...
                        'color',handles.dataaux.class_colors(class_i,:))
                end
            case 3
                if get(handles.Figures.Clusters.changevisualizationBUTTON,'value')
                    norm = zscore(handles.data.clustering_space{handles.chid});
                    GMM_PlotPrctiles(norm(plotspikes,:)',...
                        handles.dataaux.class_colors(class_i,:))
                else
                    if (size(handles.data.model{handles.chid}.mu,2)>1)
                        norm = zscore(handles.data.clustering_space{handles.chid});
                        
                        plot(norm(auxplot,:)',...
                            'color',handles.dataaux.class_colors(class_i,:))
                    else
                        plot(ones(length(auxplot),1),handles.data.clustering_space{handles.chid}(auxplot,1)',...
                            '.','color',handles.dataaux.class_colors(class_i,:))
                    end    
                end
        end

    end
end
if nclasses == 0 & handles.Figures.Clusters.mode_status==2
    plotspikes = isnan(handles.data.class_id{handles.chid});
    auxplot = find(plotspikes);
    wave_max = 1000;
    nwave = length(auxplot);
    if handles.Config.Plot==1
        if nwave>wave_max
            [v, irand]=sort(rand(1,length(auxplot)));
            auxplot = auxplot(irand(1:wave_max));
            
        end
    end
    
        
    if get(handles.Figures.Clusters.changevisualizationBUTTON,'value')
        GMM_PlotPrctiles(handles.data.waveforms{handles.chid}(plotspikes,:)',...
            [.8 .8 .8])
    else
        plot(handles.data.waveforms{handles.chid}(auxplot,:)',...
            'color',[.8 .8 .8])
    end
end
valid = ~isnan(handles.data.class_id{handles.chid});

switch handles.Figures.Clusters.mode_status
    
    case 1
        if sum(~isnan(handles.data.clustering_space{handles.chid}(valid,dimx)))>0
            if dimx~=dimy
                xlim(minmax(handles.data.clustering_space{handles.chid}(valid,dimx)'))
                ylim(minmax(handles.data.clustering_space{handles.chid}(valid,dimy)'))
            else
                xlim(minmax(handles.data.clustering_space{handles.chid}(valid,dimx)'))
                axis tight
            end
        end
        
    case 2
        if sum(valid)==0
            valid = ~valid;
        end
            ymax = max(max(handles.data.waveforms{handles.chid}(valid,:)));
            ymin = min(min(handles.data.waveforms{handles.chid}(valid,:)));
            ylim([ymin ymax])
            xlim([1 size(handles.data.waveforms{handles.chid},2)])
    case 3
        %         ymax = max(max(handles.data.clustering_space{handles.chid}(valid,:)));
        %         ymin = min(min(handles.data.clustering_space{handles.chid}(valid,:)));
        %         ylim([ymin ymax])
        
        ylim(minmax(norm(:)'))
        xlim([1 size(handles.data.clustering_space{handles.chid},2)])
        
end
box off

% axis tight

end

function [] = GMM_changestability(A,B)

GMM_plotstability(3);


end

function [] = GMM_plotstability(subplotaux)

global handles

GMM_CheckActiveElements

axes(handles.Figures.Clusters.subplots(4))
cla
axes(handles.Figures.Clusters.subplots(subplotaux))
cla

dimx = handles.Figures.Clusters.currentaxes(2);
dimy = handles.Figures.Clusters.currentaxes(1);

class_i = get(handles.Figures.Clusters.stability,'value');

% 
class_id = handles.data.class_id{handles.chid};
class_id = class_id(~isnan(class_id));
classlabels = unique(class_id);
nclasses = length(classlabels);

if(class_i<=nclasses)
    plotspikes = handles.data.class_id{handles.chid}==class_i;
    
    
%     peaktovalley = diff(minmax(handles.data.waveforms{handles.chid}(plotspikes,:))');
    peaktovalley = abs(diff(minmax(handles.data.waveforms{handles.chid}(:,:))'));
    max_diff = max(peaktovalley);
    peaktovalley = peaktovalley(plotspikes);
    
    spks = handles.data.spiketimes{handles.chid}(plotspikes)*1000;

    beginend = minmax([spks(:)]');
    timebins = 0:1000:(beginend(end));
    
    hold on
    plot(spks/60000,peaktovalley,'.',...
        'color',handles.dataaux.class_colors(class_i,:))
    hold off
%     imagesc(zscore(handles.data.clustering_space{handles.chid}(plotspikes,:))');
%     axis xy
    box off
    
    title('Cluster stability','fontsize',18)
    
    ylim([0 1.2*max_diff])
    xlim([0 max(timebins/60000)])
    ylabel('Peak-to-valley','fontsize',16)
    xlabel('Time (min)','fontsize',16)
    
    ax2=handles.Figures.Clusters.subplots(4);
    axes(ax2)

    binnedspikes = hist(spks,timebins);
    firing_rate = conv(binnedspikes,rectwin(60),'same')/60;
    hold on
    plot(timebins/60000, firing_rate, 'Parent', ax2,...
        'color',[1 1 1]*.7,...
        'linewidth',2);
    hold off
    
    box off
    set(ax2,'YColor',[1 1 1]*.7)
    ylabel('Firing rate (Hz)','fontsize',16)
    xlim([0 max(timebins/60000)])
    ylim([0 max([10 max(firing_rate)])])
    
end

end


function [ ] = GMM_PlotPrctiles( auxdata,c )

global handles


auxplotprc = prctile(auxdata',[5,25,50,75,95]);
vecwidth = [1 2 4 2 1];

for i = 1:size(auxplotprc,1)
    hold on
    plot(auxplotprc(i,:),'color',...
        c,'linewidth',...
        vecwidth(i))
    
end
hold off
axis tight
ylim(handles.Figures.Waveforms.ylim)
box off

end


function [] = GMM_startpolygon(A,B)

global handles


GMM_setselectedspks(zeros(length(handles.data.class_id{handles.chid}),1)==1);

set(handles.Figures.Clusters.startpolygon,'Enable','off') 


axes(handles.Figures.Clusters.subplots(1))
try
    [vx,vy]=getline(handles.Figures.Clusters.subplots(1),'close');
catch
    set(handles.Figures.Clusters.startpolygon,'Enable','on')
    return
end


switch handles.Figures.Clusters.mode_status
    case 1
        
        getspikesfrom_aux = get(handles.Figures.Clusters.getspikesfrom,'value');
        getspikesfrom_str = get(handles.Figures.Clusters.getspikesfrom,'string');
        getspikesfrom = getspikesfrom_str{getspikesfrom_aux};
        if strcmp(getspikesfrom,'all')
            spikes = ones(length(handles.data.class_id{handles.chid}),1);
            spikes = spikes & ~isnan(handles.data.class_id{handles.chid});
        else
            class_i = str2double(getspikesfrom);
            spikes = handles.data.class_id{handles.chid}==class_i;
        end
        
        dimx = get(handles.Figures.Clusters.display_x(1),'value');
        dimy = get(handles.Figures.Clusters.display_y(1),'value');
        
        
        xvalues = handles.data.clustering_space{handles.chid}(spikes,dimx);
        yvalues = handles.data.clustering_space{handles.chid}(spikes,dimy);
        idx = find(spikes);
        INsub = inpolygon(xvalues,yvalues,vx,vy);
        idx = idx(INsub);
        IN = spikes*0==1;
        IN(idx)=1==1;
%         plot(handles.data.clustering_space{handles.chid}(plotspikes,dimx),...
%             handles.data.clustering_space{handles.chid}(plotspikes,dimy),'.',...
%             'color',handles.dataaux.class_colors(class_i,:))

        handles.Figures.Clusters.inpolygonhandle = ...
            plot(xvalues(INsub),yvalues(INsub),'.','color',[1 1 1]*.9);
        
    case 2
        
        spikes =zeros(size(handles.data.class_id{handles.chid}));
        getspikesfrom = handles.Figures.Clusters.viewspikesVector;
        
        for class_i = getspikesfrom(:)';
            spikes = spikes | handles.data.class_id{handles.chid}==class_i;
        end
        if isempty(getspikesfrom)
            spikes = ~spikes;
        end
        %%
        yvalues = handles.data.waveforms{handles.chid}(spikes,:);
        xvalues = repmat(1:size(yvalues,2),size(yvalues,1),1);
        
        IN = linexpoly(vx,vy,xvalues,yvalues);        
        IN = find(IN==1);

        %%
        handles.Figures.Clusters.inpolygonhandle = ...
            plot(xvalues(IN,:)',yvalues(IN,:)',...
            'color',[1 1 1]*0,...
            'linewidth',2);
        
        ispikes=find(spikes==1);
        IN = ispikes(IN);
        
%         plot(handles.data.waveforms{handles.chid}(plotspikes,:)',...
%             'color',handles.dataaux.class_colors(class_i,:))
    case 3
        
        spikes =zeros(size(handles.data.class_id{handles.chid}));
        getspikesfrom = handles.Figures.Clusters.viewspikesVector;
        for class_i = getspikesfrom(:)'
            spikes = spikes | handles.data.class_id{handles.chid}==class_i;
        end
        
        
        norm_clustering_space = zscore(handles.data.clustering_space{handles.chid});        
        yvalues = norm_clustering_space(spikes,:);
        
%         yvalues = handles.data.clustering_space{handles.chid}(spikes,:);
        xvalues = repmat(1:size(yvalues,2),size(yvalues,1),1);
        IN=zeros(1,size(xvalues,1));
        for i=1:size(xvalues,1)
            if ~isempty(polyxpoly(xvalues(i,:),yvalues(i,:),vx,vy))
                IN(i) = 1;
            end
        end
        IN = find(IN==1);
        
        
        handles.Figures.Clusters.inpolygonhandle = ...
            plot(xvalues(IN,:)',yvalues(IN,:)',...
            'color',[1 1 1]*0,...
            'linewidth',2);
        
        ispikes=find(spikes==1);
        IN = ispikes(IN);
         
%         plot(handles.data.clustering_space{handles.chid}(plotspikes,:)',...
%             'color',handles.dataaux.class_colors(class_i,:))
end

% dimx = get(handles.Figures.Clusters.display_x(1),'value');
% dimy = get(handles.Figures.Clusters.display_y(1),'value');
% 
% xvalues = handles.data.clustering_space{handles.chid}(:,dimx);
% yvalues = handles.data.clustering_space{handles.chid}(:,dimy);
% 
% IN = inpolygon(xvalues,yvalues,vx,vy);
% IN = spikes & IN;

GMM_setselectedspks(IN)
% GMM_change_bck_template(handles.Figures.Clusters.selectedBKGNtemplateH)

hold(handles.Figures.Clusters.subplots(1),'on');
handles.Figures.Clusters.polygonhandle = ...
    plot(handles.Figures.Clusters.subplots(1),vx,vy,'color',[1 1 1]*.9,'linewidth',3);




% axes(handles.Figures.Clusters.subplots(2)),cla
% plot(handles.data.waveforms{handles.chid}(IN,:)','k')
% box off
% axis tight

% choosedialog

set(handles.Figures.Clusters.startpolygon,'Enable','on') 

end

function [] = GMM_setselectedspks(spks)
    global handles
    handles.Figures.Clusters.SelectedSpks = spks;
    
    try
        delete(handles.Figures.Clusters.polygonhandle)
        delete(handles.Figures.Clusters.inpolygonhandle)
    catch
        
    end

    
    GMM_change_bck_template(handles.Figures.Clusters.selectedBKGNtemplateH)
end

function [] = GMM_MoveTo(A,B)

global handles

class_id = handles.data.class_id{handles.chid};
class_id = class_id(~isnan(class_id));
classlabels = unique(class_id);

IN = handles.Figures.Clusters.SelectedSpks;
selectedoption = get(handles.Figures.Clusters.selectedBKGNtemplateH,'value');
alloptions = get(handles.Figures.Clusters.selectedBKGNtemplateH,'string');
willchange = unique(handles.data.class_id{handles.chid}(IN));

switch alloptions{selectedoption}
    case 'unsorted'
        handles.data.class_id{handles.chid}(IN) = nan;
    case 'new cluster'
        handles.data.class_id{handles.chid}(IN) = ...
            length(classlabels)+1;
        willchange=unique([willchange; length(classlabels)+1]);
    otherwise
        class = str2double(alloptions{selectedoption});
        handles.data.class_id{handles.chid}(IN) = ...
            class;
        willchange=unique([willchange; class]);
end


class_id = handles.data.class_id{handles.chid};
class_id = class_id(~isnan(class_id));
classlabels_old=classlabels;
classlabels = unique(class_id);

deleted_class=setdiff(classlabels_old,classlabels);
mode=1;
if ~isempty(deleted_class)
    handles.data.model{handles.chid}.mu=handles.data.model{handles.chid}.mu...
        (handles.data.model{handles.chid}.class~=deleted_class,:);
    handles.data.model{handles.chid}.S=handles.data.model{handles.chid}.S...
        (handles.data.model{handles.chid}.class~=deleted_class,:,:);
    handles.data.model{handles.chid}.alpha=handles.data.model{handles.chid}.alpha...
        (handles.data.model{handles.chid}.class~=deleted_class);
    
    handles.data.model{handles.chid}.class=handles.data.model{handles.chid}.class...
        (handles.data.model{handles.chid}.class~=deleted_class);
    
    mode = 0;


    clusterid = handles.data.class_id{handles.chid};
    for i = 1:length(classlabels);
        clusspikes = clusterid==classlabels(i);
        handles.data.class_id{handles.chid}(clusspikes) = i;

        handles.data.model{handles.chid}.class(handles.data.model{handles.chid}.class==classlabels(i)) = i;
    end
end
GMM_plotwaveforms(mode,willchange)
GMM_showclusters

end

function GMM_change_bck_template(A,B)

global handles


cla(handles.Figures.Clusters.subplots(2))

    
IN = handles.Figures.Clusters.SelectedSpks;
if sum(IN)>1
    
    plot(handles.Figures.Clusters.subplots(2),handles.data.waveforms{handles.chid}(IN,:)','k')
%     set(handles.Figures.Clusters.subplots(2),'box','off')
    axis(handles.Figures.Clusters.subplots(2),'tight')
    hold(handles.Figures.Clusters.subplots(2),'off')
end

optionselected = get(A,'value');
optionsstr = get(A,'string');

switch optionsstr{optionselected}
    case 'unsorted'
        
    case 'new cluster'
        
    otherwise
        class_i = str2double(optionsstr(optionselected));
    
        auxplot = handles.data.waveforms{handles.chid}';
        auxplot = ...
            auxplot(:,handles.data.class_id{handles.chid}==class_i);

        auxplotprc = prctile(auxplot',[5,25,50,75,95]);
        vecwidth = [1 2 4 2 1];

        hold(handles.Figures.Clusters.subplots(2),'on')
        for i = 1:size(auxplotprc,1)
            plot(handles.Figures.Clusters.subplots(2),auxplotprc(i,:),'color',...
                handles.dataaux.class_colors(class_i,:),'linewidth',...
                vecwidth(i))
        end
        
end

end


function GMM_runxcorr(A,B)

global handles

n1 = get(handles.Figures.Clusters.crosscorrelationN1,'value');
n2 = get(handles.Figures.Clusters.crosscorrelationN2,'value');


spiketimes = handles.data.spiketimes{handles.chid}*1000;

n1spks = spiketimes(handles.data.class_id{handles.chid}==n1);
n2spks = spiketimes(handles.data.class_id{handles.chid}==n2);

beginend = minmax([n2spks(:); n1spks(:)]');
timebins = (beginend(1)-1000):(beginend(end)+1000);

binnedspikes1 = hist(n1spks,timebins);
binnedspikes2 = hist(n2spks,timebins);

[crosscor,lags] = xcorr(binnedspikes2,binnedspikes1,30,'coeff');


axes(handles.Figures.Clusters.crosscorrelation)
cla
if n1==n2
   crosscor(lags==0)=0;
end

bar(lags,crosscor,'facecolor',handles.dataaux.class_colors(n2,:))
box off, axis tight
hold on
plot([0 0],ylim,'--','color',handles.dataaux.class_colors(n1,:),'linewidth',4)
xlabel('Time (ms)','fontsize',16)
ylabel('r','fontsize',16)
hold off

if n1==n2
   title('Auto-correlogram','fontsize',18)
else
   title(['Cross-correlogram ' int2str(n1) '-' int2str(n2)],'fontsize',18)
end

end