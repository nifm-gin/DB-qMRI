function w = abs(w)
    w = tospace(w); % only makes sens in space
    for i = 1:length(w.coarse)
        w.coarse{i} = abs(w.coarse{i});
    end
    for i = 1:length(w.wav)
        for m = 1:length(w.wav{i})
            w.wav{i}{m} = abs(w.wav{i}{m});
        end
    end
end