function GMM_savedata(A,B)

global handles

data = handles.data;

if handles.dataaux.nelectrodes>1
    nelec=handles.dataaux.nelectrodes;
    for ch_i = 1:handles.dataaux.nchannels
        data.waveforms{ch_i} = reshape(data.waveforms{ch_i},size(data.waveforms{ch_i},1),[],nelec);
    end
end

[~, fname, ext]=fileparts(handles.filename);

uisave('data',[fname '_sorting' ext])

end