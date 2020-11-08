function wz = times(wx,wy) % wx.*wy
    if any(wx.size~=wy.size)||any(wx.depth~=wy.depth)
        error('incompatible sizes: wavelet coefficients cannot be multiplied!');
    else
        if wx.fourier % multiplication makes sense for the spatial coeffs
            wx = tospace(wx);
        end
        if wy.fourier
            wy = tospace(wy);
        end
        wz = wavelet(wx);
        for i = 1:length(wz.coarse)
            wz.coarse{i} = wz.coarse{i}.*wy.coarse{i};
        end
        for i = 1:length(wz.wav)
            for m = 1:length(wz.wav{i})
                wz.wav{i}{m} = wz.wav{i}{m}.*wy.wav{i}{m};
            end
        end
    end
end