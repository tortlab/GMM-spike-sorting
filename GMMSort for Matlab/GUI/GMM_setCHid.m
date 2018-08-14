function GMM_setCHid(A,B)

global handles

chid = get(handles.Figures.Waveforms.ch_id,'value');

if isnan(chid)
    set(handles.Figures.Waveforms.ch_id,...
        'string',num2str(handles.chid));
    set(handles.Figures.Waveforms.maintext,'String',...
        'Channel id must be a number.')    
elseif ~mod(chid,1)==0
    set(handles.Figures.Waveforms.ch_id,...
        'string',num2str(handles.chid));
    set(handles.Figures.Waveforms.maintext,'String',...
        'Channel id must be an integer.')
elseif chid<=0
    set(handles.Figures.Waveforms.ch_id,...
        'string',num2str(handles.chid));
    set(handles.Figures.Waveforms.maintext,'String',...
        'Channel id must be positive.')
elseif chid>handles.dataaux.nchannels
    set(handles.Figures.Waveforms.ch_id,...
        'string',num2str(handles.chid));
    set(handles.Figures.Waveforms.maintext,'String',...
       ['We have only ' num2str(handles.dataaux.nchannels) ...
        'channels in the loaded file.'])    
else
    handles.chid = chid;
    
    handles.Figures.Clusters.viewspikesVector=unique(handles.data.class_id{handles.chid}(~isnan(handles.data.class_id{handles.chid})));

    GMM_plotwaveforms
    if get(handles.Figures.Waveforms.DispClustersTOGGLE,'value') 
        GMM_showclusters
    end
    figure(handles.Figures.Waveforms.main)
    
    set(handles.Figures.Waveforms.maintext,'String',...
        ['Displaying channel: ' ...
        num2str(handles.chid)])

end