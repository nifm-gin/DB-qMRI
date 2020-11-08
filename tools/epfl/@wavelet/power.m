function w = power(w,p) % w.^p
    w = tospace(w); % only makes sens in space
    for i = 1:length(w.coarse)
        w.coarse{i} = w.coarse{i}.^p;
    end
    for i = 1:length(w.wav)
        for m = 1:length(w.wav{i})
            w.wav{i}{m} = w.wav{i}{m}.^p;
        end
    end
end