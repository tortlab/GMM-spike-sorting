function GMM_ClusterDisplayONOFF( A,B )

global handles

get(handles.Figures.Waveforms.DispClustersTOGGLE,'value')

switch get(handles.Figures.Waveforms.DispClustersTOGGLE,'value')
    
    case false
        try
            delete(handles.Figures.Clusters.main)
            handles.Figures.Clusters = rmfield(handles.Figures.Clusters,'main');
        catch
           display('Cluster figure already closed!') 
        end
    case true
        GMM_showclusters
        
end

end