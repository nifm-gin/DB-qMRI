function w = ProxHardThresh(w,T)

if any(T<0)
    error('the threshold value should be positive!');
end

%w = tospace(w);
J = sort(w.depth);
Jm = [0,J(1:end-1)];
M = 1 + sum(double(J(:)-Jm(:)).*(2.^(numel(J):-1:1)'-1));clear Jm;
MJ = max(J);
switch numel(T)
    case 1 % threshold all the coefficients with a unique value
        for j = 1:MJ
            for s = 1:length(w.wav{j})
                w.wav{j}{s} = ProxHard(T(1),w.wav{j}{s});
            end
        end
        w.coarse{end} = ProxHard(T(1),w.coarse{end});
    case MJ+1 % threshold the coefficients with a unique value per scale
        for j = 1:MJ
            for s = 1:length(w.wav{j})
                w.wav{j}{s} = ProxHard(T(j),w.wav{j}{s});
            end
        end
        w.coarse{end} = ProxHard(T(end),w.coarse{end});
    case M % threshold the coefficients with a unique value per subband
        counter = 1;
        for j = 1:MJ
            for s = 1:length(w.wav{j})
                w.wav{j}{s} = ProxHard(T(counter),w.wav{j}{s});
                counter =counter+1;
            end
        end
        w.coarse{end} = ProxHard(T(end),w.coarse{end});
    otherwise
        error('WAVELETClass:ProxSoftThresh',...
            'Not a supported size for the threshold');
end