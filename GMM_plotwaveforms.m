function GMM_plotwaveforms( varargin )

global handles
if length(varargin)<1
    A = 0;
else
    A = varargin{1};
end

if length(varargin)<2
    B = 0;
else
    B = varargin{2};
end
    
    
    
figure(handles.Figures.Waveforms.main)

auxplot = handles.data.waveforms{handles.chid}';
handles.Figures.Waveforms.ylim = minmax(reshape(auxplot(:,~isnan(handles.data.class_id{handles.chid})),1,[]));

axes(handles.Figures.Waveforms.unsortedspikes),cla
auxplot = auxplot(:,isnan(handles.data.class_id{handles.chid}));
plot(auxplot,'color',[1 1 1]*.9)
axis tight
ylim(handles.Figures.Waveforms.ylim)
box off

classlabels = unique(handles.data.class_id{handles.chid});
classlabels = classlabels(~isnan(classlabels));
nclasses = length(classlabels);


% if ~isfield(handles.Figures.Waveforms,'cluster') || length(handles.Figures.Waveforms.cluster)~=min([7 nclasses])
if A==0
    GMM_drawwaveforms
else
    GMM_updatewaveforms(B)
end

end

function [] = GMM_drawwaveforms (A, B)
global handles

classlabels = unique(handles.data.class_id{handles.chid});
classlabels = classlabels(~isnan(classlabels));
nclasses = length(classlabels);

handles.dataaux.class_colors = linspecer(nclasses);

popupaux = num2str(classlabels(~isnan(classlabels)));
popupaux(end+1,:) = '-';

for class_i = 1:7
%     
    if class_i<=nclasses
        GMM_PlotAllSpikes( class_i )
        set( handles.Figures.Waveforms.clusterPOPUP{class_i},'String',...
            popupaux)
        set(handles.Figures.Waveforms.clusterPOPUP{class_i},'Value',...
            class_i)
        nspikes = sum(handles.data.class_id{handles.chid}==class_i);
        set(handles.Figures.Waveforms.ch_id_txtbox(class_i), 'String',...
            [num2str(nspikes) ' spikes '])
        set(handles.Figures.Waveforms.clusterNUMB{class_i},'String',int2str(class_i))
    else
        
        cla(handles.Figures.Waveforms.cluster{class_i})
        set( handles.Figures.Waveforms.clusterPOPUP{class_i},'String',...
            popupaux)
        set(handles.Figures.Waveforms.clusterPOPUP{class_i},'Value',...
            nclasses+1)
        set(handles.Figures.Waveforms.ch_id_txtbox(class_i), 'String',...
            [' '])
        
        set(handles.Figures.Waveforms.clusterNUMB{class_i},'String',' ')
    end
        
end
end



function [] = GMM_updatewaveforms(iclass)
global handles


classlabels = unique(handles.data.class_id{handles.chid});
classlabels = classlabels(~isnan(classlabels));
nclasses = length(classlabels);

handles.dataaux.class_colors = linspecer(nclasses);

popupaux = num2str(classlabels(~isnan(classlabels)));
popupaux(end+1,:) = '-';

if iclass>=1
    for class_i=iclass'
        if class_i<=7
            cla(handles.Figures.Waveforms.cluster{class_i})
            GMM_PlotAllSpikes( class_i )
            
            axis(handles.Figures.Waveforms.cluster{class_i},'tight')
            ylim(handles.Figures.Waveforms.cluster{class_i},handles.Figures.Waveforms.ylim)
            box(handles.Figures.Waveforms.cluster{class_i},'off')
            
            set(handles.Figures.Waveforms.clusterPOPUP{class_i},...
                'String', popupaux,...
                'Value',class_i)
            
            nspikes = sum(handles.data.class_id{handles.chid}==class_i);
            
            set(handles.Figures.Waveforms.ch_id_txtbox(class_i),...
                'String', [num2str(nspikes) ' spikes'])
        end
    end
    for class_i=setdiff(1:7, iclass)
        
        set( handles.Figures.Waveforms.clusterPOPUP{class_i},'String',...
            popupaux)
        set(handles.Figures.Waveforms.clusterPOPUP{class_i},'Value',...
            nclasses+1)

        set(handles.Figures.Waveforms.clusterNUMB{class_i},'String',' ')
        if class_i<=nclasses
            set(handles.Figures.Waveforms.clusterNUMB{class_i},'String',int2str(class_i))
            set(handles.Figures.Waveforms.clusterPOPUP{class_i},'Value',...
            class_i)
            set(get(handles.Figures.Waveforms.cluster{class_i},'children'),...
                'color',handles.dataaux.class_colors(class_i,:))
            nspikes = sum(handles.data.class_id{handles.chid}==class_i);
            set(handles.Figures.Waveforms.ch_id_txtbox(class_i), 'String',...
                [num2str(nspikes) ' spikes '])
        end
    end
elseif iclass~=-1

    for class_i = 1:nclasses
        if class_i<8

            if  1==1
                cla(handles.Figures.Waveforms.cluster{class_i})
                GMM_PlotAllSpikes( class_i )

                axis(handles.Figures.Waveforms.cluster{class_i},'tight')
                ylim(handles.Figures.Waveforms.cluster{class_i},handles.Figures.Waveforms.ylim)
                box(handles.Figures.Waveforms.cluster{class_i},'off')

                set(handles.Figures.Waveforms.clusterPOPUP{class_i},...
                    'String', popupaux,...
                    'Value',class_i)

                nspikes = sum(handles.data.class_id{handles.chid}==class_i);

                set(handles.Figures.Waveforms.ch_id_txtbox(class_i),...
                    'String', [num2str(nspikes) ' spikes'])
            end
        end

    end

end
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

plot(handles.Figures.Waveforms.cluster{class_i},auxplot,'color',...
    handles.dataaux.class_colors(class_i,:))
axis(handles.Figures.Waveforms.cluster{class_i},'tight');
set(handles.Figures.Waveforms.cluster{class_i},'box', 'off', 'ylim',handles.Figures.Waveforms.ylim)


end








