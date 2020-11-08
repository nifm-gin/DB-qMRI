function sref = subsref(w, S)
% Supported calls:
%
% w.field (Builtin function)
%
% w(:)          concatenation of all the coefficients in a vector
%                   (useful for norm and scalar products)
% w{j}          where j scans all subbands from fine scale to coarse
% w{j,m}        where j scans the scale and m the subbands
% w{j}{m}       idem


if numel(w) > 1
    error('WAVELETClass:subsref',...
        'Object must be scalar');
end
depth = length(S);

% detect the call of 'w(:)'
ToVector = (depth==1&&strcmp(S(1).type,'()')&&numel(S(1).subs)==1&&(strcmp(S(1).subs(1),':')));

if ToVector % provide the wavelet coefficients as a 1D vector
    undecimated_coarse = all(size(w.coarse{1})==size(w.coarse{end}));
    undecimated_wav = all(size(w.wav{end}{1})==size(w.coarse{end-1}));
    if undecimated_coarse||undecimated_wav
        ncoeff = numel(w.coarse{end});
        for i = 1:length(w.wav)
            for s = 1:length(w.wav{i})
                ncoeff = ncoeff+numel(w.wav{i}{s});
            end
        end
        sref = wavtovect(w.coarse{end},w.wav,ncoeff);
        %J = sort(w.depth);
        %Jm = [0,J(1:end-1)];%Jm = circshift(J(:),1);Jm(1) = 0;
        %M = 1 + sum(double(J(:)-Jm(:)).*(2.^(numel(J):-1:1)'-1))
        %sref = wavtovect(w.coarse{end},w.wav,M*prod(w.size));
    else
        sref = wavtovect(w.coarse{end},w.wav,prod(w.size));
    end
    return;
end

switch S(1).type
    % Use the built-in subsref for dot notation
    case '.'
        %w = tospace(w); % convert the internal representation of coefficients into space
        sref = builtin('subsref',w,S);
    case '{}'
        %w = tospace(w); % convert the internal representation of coefficients into space
        if length(S(1).subs)==1 && length(S)==1
            J = sort(w.depth);
            Jm = [0,J(1:end-1)];%Jm = circshift(J(:),1);Jm(1) = 0;
            M = 1 + sum(double(J(:)-Jm(:)).*(2.^(numel(J):-1:1)'-1));
            k = S(1).subs{1};
            if k<1||k>M
                error('WAVELETClass:subsref',...
                    'Not a supported subscripted reference');
            end
            if k==M
                sref = w.coarse{end};return;
            else
                counter = 1;
                for i = 1:length(w.wav)
                    for m = 1:length(w.wav{i})
                        if counter==k
                            sref = w.wav{i}{m};return;
                        end
                        counter = counter+1;
                    end
                end
            end
        elseif S(1).subs{1} == max(w.depth)+1,
            sref = w.coarse{end};return;
        else
            sref = w.wav{S(1).subs{1}};return;
        end
        if length(S(1).subs)>1
            sref = sref{S(1).subs{2:end}};
        end
        if length(S)>1
            sref = subsref(sref,S(2:end));
        end
    otherwise
        error('WAVELETClass:subsref',...
            'Not a supported subscripted reference');
end