function GMM_resetSorting(A,B)

global handles

handles.data.class_id{handles.chid} = ...
    handles.dataaux.class_id{handles.chid};
handles.data.model{handles.chid} = ...
    handles.dataaux.model{handles.chid};

handles.Figures.Clusters.viewspikesVector=unique(handles.data.class_id{handles.chid});

GMM_plotwaveforms
if get(handles.Figures.Waveforms.DispClustersTOGGLE,'value')
    GMM_showclusters
end
figure(handles.Figures.Waveforms.main)

end