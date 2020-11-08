function wx = plus(wx,wy)

if isa(wx, 'wavelet')&&isa(wy, 'double') % add a scalar to wavelet
    wx = tospace(wx);
    wx.coarse{end} = wx.coarse{end}+wy;
    for j = 1:length(wx.wav)
        for s =1:length(wx.wav{j})
            wx.wav{j}{s} = wx.wav{j}{s} + wy;
        end
    end
elseif isa(wx, 'double')&&isa(wy, 'wavelet') % add a scalar to wavelet
    wx = plus(wy,wx);
elseif isa(wx, 'wavelet')||isa(wy, 'wavelet') % add two wavelet objects
    if any(wx.size~=wy.size)||any(wx.depth~=wy.depth)
        error('incompatible parameters: wavelet coefficients cannot be added!');
    end
    if wy.fourier~=wx.fourier
        wx = tospace(wx);
        wy = tospace(wy);
    end
    wx.coarse{end} = wx.coarse{end}+wy.coarse{end};
    wx.wav = wavadd(wx.wav,wy.wav,1);
else
    error('inputs should be both wavelet coefficients');
end

