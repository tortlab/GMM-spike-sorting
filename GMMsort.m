function [varargout] = GMMsort(varargin)
% GMMsort Implements spike sorting using Gaussian mixture models (GMMs)
% 
%     [] = GMMsort opens the GUI to perform the classification
% 
%     [MODEL, CLUSTERID, FEATURES] = GMMsort(X), where X is a matrix
%     containing the waveforms to be sorted (each line is a waveform),
%     returns the fitted GMM into MODEL and the cluster identity of each
%     waveform in the vector CLUSTERID.
% 
%     [MODEL, CLUSTERID, FEATURES] = GMMsort(X,'nGauss', NGAUSS) where
%     NGAUSS is the number of Gaussians used to compute the unidimensional
%     GMMs.
% 
%     [MODEL, CLUSTERID, FEATURES] = GMMsort(X,'nGaussClus', NGAUSSCLUS)
%     where NGAUSSCLUS is the number of Gaussians used to compute the
%     multidimensional GMMs, and thus, the maximum number of possible
%     clusters.
% 
%     [MODEL, CLUSTERID, FEATURES] = GMMsort(X,'replicates', REPLICATES)
%     where REPLICATES is the number of times each unimodal GMM is repeated
%     to estimate the cluster separability of the feature.
% 
%     [] = GMMsort('Option', OPTION) opens the
%     GUI with the specified option value ('nGauss','nGaussClus',
%     'replicates')
% 
% Example of GUI usage: Run 'GMMsort'. 'Load' the file sample_waveforms.mat
% and 'Run Sorting'. Edit the classification in the 'Clusters' window.
% 
% 'sample_waveforms.mat' contains the spike times and waveforms of the
% reccorded neurons in a variable struct called 'data' with the following
% fields:
%   data.waveforms - a 1 by NCH cell vector, where NCH is the number of
%       channels. Each cell must have a NSPK by N matrix, with NSPK
%       waveforms of length N.
% 
%   data.spiketimes - a 1 by NCH cell vector, where NCH is the number of
%       channels. Each cell must have a NSPK by 1 matrix, with NSPK spike
%       times.
% 
% See detailed explanation of the GUI in the instruction file:
% 'GMMinstruction.pdf'.
% 
% B. C. Souza January, 2018
% Brain Institute, Natal, Brazil



%% SET PARAMETERS

global parameters
parameters =[];

% maximum number of gaussians to fit in one dimension.
parameters.maxGauss         = 8;

% maximum number of iterations of (exp maximization) EM algorithm
parameters.optgmfit.max_iter    = 10000;
parameters.optgmfit.conv_factor = 1e-6;

% number of models to calculate to check robustness
parameters.nof_replicates	= 10;

% maximum number of gaussians to overfit in multidim space. 12 for
% single-wire. 20 for tetrodes.
parameters.ngaussovfit      = 12;


if ~isempty(varargin)
    [REG,prop]=parseparams(varargin(2:end));
    
    idx=find(strcmpi('nGauss',prop));
    
    if ~isempty(idx)
        parameters.maxGauss=prop{idx+1};
    end
    
    idx=find(strcmpi('nGaussClus',prop));
    if ~isempty(idx)
        parameters.ngaussovfit=prop{idx+1};
    end
    
    idx=find(strcmpi('replicates',prop));
    if ~isempty(idx)
        parameters.nof_replicates=prop{idx+1};
    end
end
%%
varargout={};
if mod(length(varargin),2)==0
    GMMsort_gui
else
    waveforms = varargin{1};
    
    features = extract_features(waveforms,parameters);
    
    [model,clusterid] = clusterize(features(:,1:5),parameters);
    
    varargout{1} = model;
    varargout{2} = clusterid;
    varargout{3} = features;
    
    % plot the model elipses and features in the first 2 dimensions..
    dim= [1 2];
    plot_model(model, clusterid, features, dim)
end

end

function GMMsort_gui()
%% creates figure for waveforms
global handles
handles =[];

% % 
% % handles.newdata = 0;
handles.data.waveforms = {};
% % handles.data.class_id = {};
% % 
handles.chid = 0;


handles.Figures.Waveforms.main = figure('units','normalized','position',[.1 .1 .7 .8],'DeleteFcn',@GMM_close);clf
set(gcf,'color','w')
set(handles.Figures.Waveforms.main, 'MenuBar', 'none');

handles.Figures.Waveforms.LoadData = ...
    uicontrol ('Style', 'pushbutton', 'String', 'LOAD','Units','normalized',...
    'Position', ...
    [0.01 0.94 0.05 0.04],...
    'Callback', @GMM_loaddata,...
    'BackgroundColor', [1 1 1]*.9,'ForegroundColor', [1 1 1]*0);
set(handles.Figures.Waveforms.LoadData,'FontName','Arial')
set(handles.Figures.Waveforms.LoadData,'fontsize',17)

handles.Figures.Waveforms.SaveData = ...
    uicontrol ('Style', 'pushbutton', 'String', 'SAVE','Units','normalized',...
    'Position', ...
    [0.065 0.94 0.05 0.04],...
    'Callback', @GMM_savedata,...
    'BackgroundColor', [1 1 1]*.9,'ForegroundColor', [1 1 1]*0,...
    'FontName','Arial',...
    'fontsize',17);

maintextbox = [0.17 .93 .8 .04];
handles.Figures.Waveforms.maintext = ...
    uicontrol('Style','text','Units','normalized',...
    'Position',...
    maintextbox,...
    'String','Please load some data... ',...
    'BackgroundColor', [1 1 1]*.99,...
    'ForegroundColor', [1 1 1]*0,...
    'horizontalAlignment', 'left',...
    'fontsize',15);

positionaux = [0.8 .93 .19 .04];
handles.Figures.Waveforms.filename = ...
    uicontrol('Style','text','Units','normalized',...
    'Position',...
    positionaux,...
    'String',' No file loaded. ',...
    'BackgroundColor', [1 1 1]*.99,...
    'ForegroundColor', [1 1 1]*0,...
    'horizontalAlignment', 'right',...
    'fontsize',15);

handles.Figures.Waveforms.ch_id_txtbox = uicontrol('Style', 'text',...
    'Units','normalized',...
    'String','Ch id: ',...
    'Position', [.01 .9 .08 .02],...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'horizontalAlignment', 'left',...
    'fontsize',14);

handles.Figures.Waveforms.ch_id = ...
            uicontrol('Style', 'popup',...
            'Units','normalized',...
            'String', '-',...
            'Position', [.05 .895 .025 .0275],...
            'BackgroundColor', [1 1 1]*1,...
            'ForegroundColor', [1 1 1]*0,...
            'Callback',@GMM_setCHid,...
            'horizontalAlignment', 'center',...
            'fontsize',14);

positionaux = [0.0825 0.8875 0.08 0.04];
handles.Figures.Waveforms.runSorting = uicontrol ('Style', 'pushbutton', ...
    'String', 'Run Sorting',...
    'Units','normalized',...
    'Position', ...
    [0.0825 0.8875 0.08 0.04],...
    'Callback', @GMM_runSorting,...
    'BackgroundColor', [1 1 1]*.9,'ForegroundColor', [1 1 1]*0,...
    'fontsize',14);

positionaux(1) = sum(positionaux([1,3]))+0.01;
handles.Figures.Waveforms.resetSorting = uicontrol ('Style', 'pushbutton', ...
    'String', 'Reset',...
    'Units','normalized',...
    'Position', ...
    positionaux,...
    'Callback', @GMM_resetSorting,...
    'BackgroundColor', [1 1 1]*.9,'ForegroundColor', [1 1 1]*0,...
    'fontsize',14);

positionaux(1) = sum(positionaux([1,3]))+0.01;
handles.Figures.Waveforms.DispClustersTOGGLE = ...
    uicontrol('Style', 'togglebutton',...
    'Units','normalized',...
    'Position', positionaux,...
    'String','Clusters',...
    'BackgroundColor', [1 1 1]*.9,...
    'ForegroundColor', [1 1 1]*0,...
    'Callback',@GMM_ClusterDisplayONOFF,...
    'fontsize',14,...
    'horizontalAlignment', 'center');

positionaux(1) = sum(positionaux([1,3]))+0.01;
handles.Figures.Waveforms.RunSortingAllchs = ...
    uicontrol('Style', 'pushbutton',...
    'Units','normalized',...
    'Position', positionaux,...
    'String','Run ALL Chs',...
    'BackgroundColor', [1 1 1]*.9,...
    'ForegroundColor', [1 1 1]*0,...
    'Callback',@GMM_SSallCHs,...
    'fontsize',14,...
    'horizontalAlignment', 'center');

positionaux(1) = sum(positionaux([1,3]))+0.01;
% positionaux = [0.245 0.8875 0.07 0.04];
handles.Figures.Waveforms.EditUnsorted = ...
    uicontrol('Style', 'pushbutton',...
    'Units','normalized',...
    'Position', positionaux,...
    'String','Edit Unsorted',...
    'BackgroundColor', [1 1 1]*.9,...
    'ForegroundColor', [1 1 1]*0,...
    'Callback',@GMM_EditUnsorted,...
    'fontsize',14,...
    'horizontalAlignment', 'center');

positionaux(1) = sum(positionaux([1,3]))+0.01;
% positionaux = [0.245 0.8875 0.07 0.04];
handles.Figures.Waveforms.DiscardUnsorted = ...
    uicontrol('Style', 'pushbutton',...
    'Units','normalized',...
    'Position', positionaux,...
    'String','Discard Unsorted',...
    'BackgroundColor', [1 1 1]*.9,...
    'ForegroundColor', [1 1 1]*0,...
    'Callback',@GMM_DiscardUnsorted,...
    'fontsize',14,...
    'horizontalAlignment', 'center');

positionaux(1) = sum(positionaux([1,3]))+0.01;
% positionaux = [0.105+class_i*0.25 .875 .02 .026];
% positionaux=positionaux-positionmod;
handles.Config.Plot = 1;
handles.Figures.Waveforms.PlotAllTOGGLE = ...
    uicontrol('Style', 'togglebutton',...
    'String','Plot all waveforms',...
    'Value',handles.Config.Plot,...
    'Units','normalized',...
    'Position', positionaux,...
    'BackgroundColor', [1 1 1]*1,...
    'ForegroundColor', [1 1 1]*0,...
    'Callback',@GMM_ChangePlotAll,...
    'horizontalAlignment', 'center');


positionaux = [0.005 .46 .22 .4];
handles.Figures.Waveforms.UnsortedPanel = ...
    uipanel('Title','Unsorted Waveforms','FontSize',13,...
    'BackgroundColor',[1 1 1]*.97,...
    'bordertype','etchedout',...
    'ShadowColor',[1 1 1]*.8,...
    'highlightcolor',[1 1 1]*.5,...
    'Position',positionaux);
%     'Position',[.0 .015 .275 .45]);

positionaux = [0.025 .5 .18 .325];
% positionaux = [0.015 .5 .2 .375];
handles.Figures.Waveforms.unsortedspikes = ...
    subplot('Position',positionaux);



%%


for class_i = 1:7
    
    if class_i<4
        positionmod=[0.03 0.03 0 0];
        
        positionaux = [0.05+class_i*0.25 .5 .22 .375];
        positionaux=positionaux-positionmod;
        handles.Figures.Waveforms.cluster{class_i} = ...
            subplot('Position',positionaux);
        
        positionaux = [0.055+class_i*0.25 .875 .025 .026];
        positionaux=positionaux-positionmod;
        handles.Figures.Waveforms.clusterNUMB{class_i} = ...
            uicontrol('Style', 'edit',...
            'Units','normalized',...
            'String', '',...
            'Position', positionaux,...
            'BackgroundColor', [1 1 1]*1,...
            'ForegroundColor', [1 1 1]*0,...
            'Callback',@GMM_ChangeClusterNumber,...
            'horizontalAlignment', 'center');
        set(handles.Figures.Waveforms.clusterNUMB{class_i},'String',...
            num2str(class_i)) 
        
        positionaux(1) = sum(positionaux([1,3])) + .005;
        handles.Figures.Waveforms.clusterPOPUP{class_i} = ...
            uicontrol('Style', 'popup',...
            'Units','normalized',...
            'String', ' ',...
            'Position', positionaux+[0 -0.015 0 .015],...
            'BackgroundColor', [1 1 1]*1,...
            'ForegroundColor', [1 1 1]*0,...
            'Callback',@GMM_ChangeDisplayedCluster,...
            'horizontalAlignment', 'center');
%         set(handles.Figures.Waveforms.clusterPOPUP{class_i},'Value',...
%             class_i) 
        
        positionaux(1) = sum(positionaux([1,3])) + .005;
        handles.Figures.Waveforms.clusterTOGGLE{class_i} = ...
            uicontrol('Style', 'togglebutton',...
            'String','%',...
            'Units','normalized',...
            'Position', positionaux,...
            'BackgroundColor', [1 1 1]*1,...
            'ForegroundColor', [1 1 1]*0,...
            'Callback',@GMM_ChangeVisualization,...
            'horizontalAlignment', 'center');
        
        positionaux(1) = sum(positionaux([1,3])) + .005;
        handles.Figures.Waveforms.delBUTTON{class_i} = ...
            uicontrol('Style', 'pushbutton',...
            'Units','normalized',...
            'Position', positionaux,...
            'String','DEL',...
            'BackgroundColor', [1 .5 .5]*.8,...
            'ForegroundColor', [1 1 1]*0,...
            'Callback',@GMM_deletecluster,...
            'horizontalAlignment', 'center');
        
        positionaux = [0.19+class_i*0.25 .8775 .07 .02];
        positionaux=positionaux-positionmod;
%         nspikes = sum(handles.data.class_id{handles.chid}==class_i);
        handles.Figures.Waveforms.ch_id_txtbox(class_i) = ...
            uicontrol('Style', 'text',...
            'Units','normalized',...
            'String',' ',...
            'Position', positionaux,...
            'BackgroundColor', [1 1 1]*1,...
            'ForegroundColor', [1 1 1]*0,...
            'horizontalAlignment', 'right',...
            'fontsize',14);
        
    elseif class_i>=4 && class_i<8
        positionmod(2)=0.015;
        
        positionaux = [0.05+(class_i-4)*0.25 .05 .22 .375];
        positionaux=positionaux-positionmod;
        handles.Figures.Waveforms.cluster{class_i} = ...
            subplot('Position',positionaux);
        
        
        positionaux = [0.055+(class_i-4)*0.25 .425 .025 .026];
        positionaux=positionaux-positionmod;
        handles.Figures.Waveforms.clusterNUMB{class_i} = ...
            uicontrol('Style', 'edit',...
            'Units','normalized',...
            'String', '',...
            'Position', positionaux,...
            'BackgroundColor', [1 1 1]*1,...
            'ForegroundColor', [1 1 1]*0,...
            'Callback',@GMM_ChangeClusterNumber,...
            'horizontalAlignment', 'center');
        set(handles.Figures.Waveforms.clusterNUMB{class_i},'String',...
            num2str(class_i)) 
        
        positionaux(1) = sum(positionaux([1,3])) + .005;
        handles.Figures.Waveforms.clusterPOPUP{class_i} = ...
            uicontrol('Style', 'popup',...
            'Units','normalized',...
            'String', ' ',...
            'Position', positionaux,...
            'BackgroundColor', [1 1 1]*1,...
            'ForegroundColor', [1 1 1]*0,...
            'Callback',@GMM_ChangeDisplayedCluster,...
            'horizontalAlignment', 'center');
        
        positionaux(1) = sum(positionaux([1,3])) + .005;
        handles.Figures.Waveforms.clusterTOGGLE{class_i} = ...
            uicontrol('Style', 'togglebutton',...
            'String','%',...
            'Units','normalized',...
            'Position', positionaux,...
            'BackgroundColor', [1 1 1]*1,...
            'ForegroundColor', [1 1 1]*0,...
            'Callback',@GMM_ChangeVisualization,...
            'horizontalAlignment', 'center');
        
        positionaux(1) = sum(positionaux([1,3])) + .005;
        handles.Figures.Waveforms.delBUTTON{class_i} = ...
            uicontrol('Style', 'pushbutton',...
            'Units','normalized',...
            'Position', positionaux,...
            'String','DEL',...
            'BackgroundColor', [1 .5 .5]*.8,...
            'ForegroundColor', [1 1 1]*0,...
            'Callback',@GMM_deletecluster,...
            'horizontalAlignment', 'center');
        
        positionaux = [0.19+(class_i-4)*0.25 .4275 .07 .02];
        positionaux=positionaux-positionmod;
        handles.Figures.Waveforms.ch_id_txtbox(class_i) = ...
            uicontrol('Style', 'text',...
            'Units','normalized',...
            'String',' ',...
            'Position', positionaux,...
            'BackgroundColor', [1 1 1]*1,...
            'ForegroundColor', [1 1 1]*0,...
            'horizontalAlignment', 'right',...
            'fontsize',14);
        
    end
end

end

function GMM_close(A,B)
global handles
    try
        close(handles.Figures.Clusters.main)
    catch e 
        
    end
clear handles
end


function GMM_SSallCHs(A,B)

global handles

for chi = 1:handles.dataaux.nchannels
    
    handles.chid = chi;
%     set(handles.Figures.Waveforms.ch_id,'String',...
    set(handles.Figures.Waveforms.ch_id,'value',...
        num2str(handles.chid))
    
    GMM_runSorting
    pause(.1)
    
end

end

function [ ] = GMM_ChangeClusterNumber( A,B )

global handles

aux = cell2mat(handles.Figures.Waveforms.clusterNUMB)==A;
new_class = str2num(get(handles.Figures.Waveforms.clusterNUMB{aux},'string'));
old_class = (get(handles.Figures.Waveforms.clusterPOPUP{aux},'value'));

new_class_idx = handles.data.class_id{handles.chid}==new_class;
old_class_idx = handles.data.class_id{handles.chid}==old_class;


handles.data.class_id{handles.chid}(old_class_idx) = new_class;
handles.data.class_id{handles.chid}(new_class_idx) = old_class;

model_old = find(handles.data.model{handles.chid}.class==old_class);
model_new = find(handles.data.model{handles.chid}.class==new_class);

handles.data.model{handles.chid}.class(model_old) = new_class;
handles.data.model{handles.chid}.class(model_new) = old_class;


% cellfun(@(x) str2num(get(x,'string')) ,handles.Figures.Waveforms.clusterNUMB)
aux_new = find(cellfun(@(x) (get(x,'value')) ,handles.Figures.Waveforms.clusterPOPUP)==new_class);
aux_old = find(cellfun(@(x) (get(x,'value')) ,handles.Figures.Waveforms.clusterPOPUP)==old_class);

cellfun(@(x) (set(x,'string',num2str(old_class))),handles.Figures.Waveforms.clusterNUMB(aux_old))
cellfun(@(x) (set(x,'string',num2str(new_class))),handles.Figures.Waveforms.clusterNUMB(aux_new))
aux = [aux_new aux_old];
for iaux=aux(:)'
    GMM_ChangeVisualization(handles.Figures.Waveforms.clusterTOGGLE{iaux})
end

axis tight
ylim(handles.Figures.Waveforms.ylim)
box off

end


function [ ] = GMM_ChangeDisplayedCluster( A,B )

global handles

aux = cell2mat(handles.Figures.Waveforms.clusterPOPUP)==A;
axes(handles.Figures.Waveforms.cluster{aux}),cla

GMM_ChangeVisualization(handles.Figures.Waveforms.clusterTOGGLE{aux})
set(handles.Figures.Waveforms.clusterNUMB{aux},'String',num2str(get(handles.Figures.Waveforms.clusterPOPUP{aux},'value')))

axis tight
ylim(handles.Figures.Waveforms.ylim)
box off

end


function [ ] = GMM_ChangeVisualization( A,B )

global handles

aux = cell2mat(handles.Figures.Waveforms.clusterTOGGLE)==A;
class_i = get(handles.Figures.Waveforms.clusterPOPUP{aux},'value');

axes(handles.Figures.Waveforms.cluster{aux})
cla
if class_i>size(handles.dataaux.class_colors,1)
    return;
else

    switch get(A,'value')

        case false
            GMM_PlotAllSpikes( class_i )
            set(A,'String','%')

        case true
            GMM_PlotPrctiles( class_i )
            set(A,'String','=')
    end
end
end

function [ ] = GMM_PlotPrctiles( class_i )

global handles

auxplot = handles.data.waveforms{handles.chid}';
auxplot = auxplot(:,handles.data.class_id{handles.chid}==class_i);

auxplotprc = prctile(auxplot',[5,25,50,75,95]);
vecwidth = [1 2 4 2 1];

for i = 1:size(auxplotprc,1)
    hold on
    plot(auxplotprc(i,:),'color',...
        handles.dataaux.class_colors(class_i,:),'linewidth',...
        vecwidth(i))
    
end
hold off
axis tight
ylim(handles.Figures.Waveforms.ylim)
box off

end


function [ ] = GMM_PlotAllSpikes( class_i )

global handles

class_idx = handles.data.class_id{handles.chid}==class_i;
auxplot = handles.data.waveforms{handles.chid}';
auxplot = auxplot(:,class_idx);

wave_max = 1000;
if handles.Config.Plot==1
    if size(auxplot,2)>wave_max
        [v, irand]=sort(rand(1,size(auxplot,2)));
        auxplot = auxplot(:,irand(1:wave_max));
    end
end

plot(auxplot,'color',...
    handles.dataaux.class_colors(class_i,:))
axis tight
ylim(handles.Figures.Waveforms.ylim)
box off

end


function [ ] = GMM_deletecluster(A,B )

global handles

aux = cell2mat(handles.Figures.Waveforms.delBUTTON)==A;
class_i = get(handles.Figures.Waveforms.clusterPOPUP{aux},'value');

handles.data.model{handles.chid}.class(handles.data.model{handles.chid}.class==class_i)=nan;

clusterid = handles.data.class_id{handles.chid};
clusterid(clusterid==class_i)=nan;

clusterlabels = unique(clusterid(~isnan(clusterid)));

for i = 1:length(clusterlabels);
    
    clusspikes = clusterid==clusterlabels(i);
    handles.data.class_id{handles.chid}(clusspikes) = i;
    handles.data.model{handles.chid}.class(handles.data.model{handles.chid}.class==clusterlabels(i))=i;
end
handles.data.class_id{handles.chid}(isnan(clusterid))=nan;

GMM_plotwaveforms(0)

if get(handles.Figures.Waveforms.DispClustersTOGGLE,'value')
    GMM_showclusters
end

end


function GMM_ChangePlotAll(A,B)
    global handles
    if handles.Config.Plot == 1
        handles.Config.Plot = 0;
        set(handles.Figures.Waveforms.PlotAllTOGGLE,'String','Plot <1000');
    else
        handles.Config.Plot = 1;
        set(handles.Figures.Waveforms.PlotAllTOGGLE,'String','Plot all');
    end
    GMM_plotwaveforms
    if get(handles.Figures.Waveforms.DispClustersTOGGLE,'value')
        GMM_showclusters
    end
end


function GMM_EditUnsorted(A,B)

global handles

clusterid = handles.data.class_id{handles.chid};
clusterid = clusterid+1;

clusterid(isnan(clusterid)) = 1;
% clusterid(clusterid==class_i)=nan;

handles.data.class_id{handles.chid}=clusterid;
if isfield(handles.data.model{handles.chid},'class')
    handles.data.model{handles.chid}.class= handles.data.model{handles.chid}.class+1;
end
GMM_plotwaveforms

if get(handles.Figures.Waveforms.DispClustersTOGGLE,'value')
    GMM_showclusters
end

end


function GMM_DiscardUnsorted(A,B)
global handles


clusterid = handles.data.class_id{handles.chid};
keep = find(~isnan(clusterid));


handles.data.waveforms{handles.chid} = handles.data.waveforms{handles.chid}(keep,:);
handles.data.class_id{handles.chid} = handles.data.class_id{handles.chid}(keep,:);
handles.dataaux.class_id{handles.chid} = handles.dataaux.class_id{handles.chid}(keep,:);
handles.data.clustering_space{handles.chid} = handles.data.clustering_space{handles.chid}(keep,:);
handles.data.classification_uncertainty{handles.chid} = handles.data.classification_uncertainty{handles.chid}(keep,:);
handles.data.spiketimes{handles.chid} = handles.data.spiketimes{handles.chid}(keep,:);

GMM_plotwaveforms(1,-1)

if get(handles.Figures.Waveforms.DispClustersTOGGLE,'value')
    GMM_showclusters
end

end



