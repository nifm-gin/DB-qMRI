function wx = minus(wx,wy)

if ~isa(wx, 'wavelet')||~isa(wy, 'wavelet') % add two wavelet objects
    error('inputs should be both wavelet coefficients');
end

if any(wx.size~=wy.size)||any(wx.depth~=wy.depth)
    error('incompatible parameters: wavelet coefficients cannot be added!');
end

if wy.fourier~=wx.fourier
    wx = tospace(wx);
    wy = tospace(wy);
end
wx.coarse{end} = wx.coarse{end}-wy.coarse{end};
wx.wav = wavadd(wx.wav,wy.wav,-1);