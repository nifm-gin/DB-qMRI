function wz = rdivide(wx,wy) % wx./wy
if isa(wx, 'wavelet')&&isa(wy, 'wavelet') % add two wavelet objects
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
        for i = 1:length(wz.coarse)
            wz.coarse{i} = wz.coarse{i}./wy.coarse{i};
        end
        for i = 1:length(wz.wav)
            for m = 1:length(wz.wav{i})
                wz.wav{i}{m} = wz.wav{i}{m}./wy.wav{i}{m};
            end
        end
    end
elseif any([isa(wx, 'double'), isa(wy, 'double')]) % add a constant to a wavelet
    if isa(wx,'double') % case: scalar./wavelet
        alpha = wx;
        wz = tospace(wy); % operation only makes sense in space for the spatial coeffs
        wz.coarse{end} = alpha./wz.coarse{end};
        for i = 1:length(wz.wav)
            for m = 1:length(wz.wav{i})
                wz.wav{i}{m} = alpha./wz.wav{i}{m};
            end
        end
    else % case: wavelet ./ scalar
        wz = (1/wy)*wx;
    end
else
    error('addition with wavelet: wrong type!')
end
end