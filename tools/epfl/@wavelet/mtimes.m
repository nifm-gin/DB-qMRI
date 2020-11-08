function w = mtimes(alpha,w) % alpha*wx
% to define the product of wavelet coeff with a scalar
if isa(alpha, 'double')&&isa(w, 'wavelet') % case scalar * wavelet
    J = sort(w.depth(:));
    Jm = [0,J(1:end-1)];% circshift(J(:),1);Jm(1) = 0;
    M = 1 + sum(double(J(:)-Jm(:)).*(2.^(numel(J):-1:1)'-1));
    MJ = J(end);
    
    if all(numel(alpha)~=[1,MJ+1,M])
        error('WAVELETClass:mtimes',...
            'Not a supported size for multiplication');
    end
    w.wav = wavmult(alpha,w.wav);
    w.coarse{end} = alpha(end)*w.coarse{end};
    %     switch numel(alpha)
    %         case 1 % multiply with a unique scalar
    %             for i = 1:MJ
    %                 for m = 1:length(w.wav{i})
    %                     w.wav{i}{m} = alpha*w.wav{i}{m};
    %                 end
    %             end
    %         case MJ+1 % multiply with a scalar per scale
    %             for i = 1:MJ
    %                 for m = 1:length(w.wav{i})
    %                     w.wav{i}{m} = alpha(i)*w.wav{i}{m};
    %                 end
    %             end
    %         case M % multiply with a scalar per subband
    %             counter = 1;
    %             for i = 1:MJ
    %                 for m = 1:length(w.wav{i})
    %                     w.wav{i}{m} = alpha(counter)*w.wav{i}{m};
    %                     counter =counter+1;
    %                 end
    %             end
    %     end
elseif isa(alpha, 'wavelet')&&isa(w, 'double') % case wavelet * scalar
    w = mtimes(w,alpha);
elseif isa(alpha, 'wavelet')&&isa(w, 'wavelet') % case wavelet * scalar
    w = times(w,alpha);
else
    error('WAVELETClass:mtimes',...
        'Not a supported class for multiplication with wavelet');
end