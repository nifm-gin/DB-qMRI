function w = uminus(w) % -w
    for i = 1:length(w.coarse)
        w.coarse{i} = -w.coarse{i};
    end
    for i = 1:length(w.wav)
        for m = 1:length(w.wav{i})
            w.wav{i}{m} = -w.wav{i}{m};
        end
    end
end