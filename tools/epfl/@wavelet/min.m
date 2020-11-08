function wz = min(wx,wy)
if isa(wx, 'wavelet')&&isa(wy, 'wavelet') % max on two wavelet objects
    if any(wx.size~=wy.size)||any(wx.depth~=wy.depth)
        error('incompatible sizes: wavelet coefficients cannot be multiplied!');
    else
        if wx.fourier % operation makes sense for the spatial coeffs
            wx = tospace(wx);
        end
        if wy.fourier
            wy = tospace(wy);
        end
        wz = wavelet(wx);
        
        for i = 1:length(wz.wav)
            for m = 1:length(wz.wav{i})
                wz.wav{i}{m} = min(wx.wav{i}{m},wy.wav{i}{m});
            end
        end
        wz.coarse{i} = min(wx.coarse{i},wy.coarse{i});
    end
elseif any([isa(wx, 'double'), isa(wy, 'double')]) % add a constant to a wavelet
    if isa(wx,'double') % case: max(scalar,wavelet)
        alpha = wx;
        wz = tospace(wy); % operation only makes sense in space for the spatial coeffs
        for i = 1:length(wz.wav)
            for m = 1:length(wz.wav{i})
                wz.wav{i}{m} = min(alpha,wy.wav{i}{m});
            end
        end
        wz.coarse{end} = min(alpha,wz.coarse{end});
    else % case: max(wavelet,scalar)
        wz = min(wy,wx);
    end
else
    error('addition with wavelet: wrong type!')
end
end