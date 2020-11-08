function w = gt(w,alpha)
% to define the comparison of wavelet coeff with a scalar
if isa(alpha, 'double')&&isa(w, 'wavelet') % case wavelet > scalar
    if w.fourier
        w = tospace(w);
    end
    w.coarse{end} = double(abs(w.coarse{end})>alpha);
    for i = 1:length(w.wav)
        for m = 1:length(w.wav{i})
            w.wav{i}{m} = double(abs(w.wav{i}{m})>alpha);
        end
    end
    %w.wav = wavgt(alpha,w.wav);
    %w.coarse{end} = (alpha<abs(w.coarse{end}));
elseif isa(alpha, 'wavelet')&&isa(w, 'double') % case  scalar > wavelet
    wz = gt(alpha,w);
else
    error('WAVELETClass:mtimes',...
        'Not a supported class for multiplication with wavelet');
end