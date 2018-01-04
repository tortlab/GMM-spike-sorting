function GMM_loaddata(A,B)

global handles

handles.newdata = 0;
try
    [filename, path]= uigetfile;
    if filename~=0
        
        load(fullfile(path,filename));
        handles.data.waveforms = data.waveforms;
        handles.dataaux.nchannels = length(data.waveforms);
        handles.newdata = 1;
        handles.filename = filename;
        
        try
            handles.data = rmfield(handles.data,'class_id');
            handles.data = rmfield(handles.data,'clustering_space');
            handles.data = rmfield(handles.data,'classification_uncertainty');
        catch e
            
        end
        
        
        if ~isfield(data,'class_id')
            handles.data.class_id = cell(1,handles.dataaux.nchannels);
            for ch_i = 1:handles.dataaux.nchannels
                handles.data.class_id{ch_i} = ...
                    nan(size(handles.data.waveforms{ch_i},1),1);
            end
            handles.dataaux.class_id = handles.data.class_id;
        else
            handles.data.class_id = data.class_id;
            handles.dataaux.class_id = handles.data.class_id;
        end
        
        if ~isfield(data,'clustering_space')
            handles.data.clustering_space = handles.data.class_id;
            handles.data.classification_uncertainty = ...
                handles.data.class_id;
        else
            handles.data.clustering_space = data.clustering_space;
            handles.data.classification_uncertainty = data.classification_uncertainty;
        end
        
%         if ~isfield(data,'maxuncertainty')
%             handles.data.maxuncertainty = 5*ones(handles.dataaux.nchannels,1);
%             handles.dataaux.maxuncertainty = handles.data.maxuncertainty;
%         else
%             handles.data.maxuncertainty = data.maxuncertainty;
%             handles.dataaux.maxuncertainty = data.maxuncertainty;
%         end
        
        if ~isfield(data,'model')
            handles.data.model = cell(1,handles.dataaux.nchannels);
            handles.dataaux.model = handles.data.model;
        else
            handles.data.model = data.model;
            handles.dataaux.model = handles.data.model;
        end
        
        if ~isfield(data,'spiketimes')
            handles.data.spiketimes = cell(1,handles.dataaux.nchannels);
            for ch_i = 1:handles.dataaux.nchannels
                handles.data.spiketimes{ch_i} = ...
                    zeros(size(handles.data.waveforms{ch_i},1),1);
            end
        else
            handles.data.spiketimes = data.spiketimes;
        end
        
        
        handles.dataaux.nelectrodes = size(handles.data.waveforms{1},3);
        if handles.dataaux.nelectrodes>1
            nelec=handles.dataaux.nelectrodes;
            for ch_i = 1:handles.dataaux.nchannels
                concatwaveforms=[];
                for ielec = 1:nelec
                    concatwaveforms = [concatwaveforms squeeze(handles.data.waveforms{ch_i}(:,:,ielec))];
                end
                handles.data.waveforms{ch_i} = concatwaveforms;
            end
        end
    end
catch errormsg
    display(errormsg.message)
end

if isempty(handles.data.waveforms)
    set(handles.Figures.Waveforms.maintext ,'String',...
        'Waveforms not loaded.')
elseif ~iscell(handles.data.waveforms)
    set(handles.Figures.Waveforms.maintext ,'String',...
        'Waveforms file should be a cell.')       
else
    if handles.newdata
        
        msgaux = [num2str(handles.dataaux.nchannels) ' channels loaded.'];
        
        handles.chid = 1;
        set(handles.Figures.Waveforms.maintext,'String',msgaux)
        set(handles.Figures.Waveforms.filename,'String',handles.filename)
        
        popupaux= cellfun(@num2str,num2cell(1:handles.dataaux.nchannels),'Uniformoutput',0);
        set(handles.Figures.Waveforms.ch_id,'String', popupaux,'Value',handles.chid)
        
        GMM_plotwaveforms
        if get(handles.Figures.Waveforms.DispClustersTOGGLE,'value')
            GMM_showclusters
        end

    end
    
end

end