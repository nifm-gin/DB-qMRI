function w = subsasgn(w, S, val)
% Supported calls:
%
% w.field = val (Builtin function)
%
% w{j} = tab        where j scans all subbands from fine scale to coarse
% w{j,m} = tab      where j scans the scale and m the subbands
% w{j}{m} = tab     idem

depth = length(S);

if numel(w) > 1
    error('WAVELETClass:subsref',...
        'Object must be scalar');
end
switch S(1).type
    % Use the built-in subsref for dot notation
    case '.'
        %w = tospace(w); % convert the internal representation of coefficients into space
        w = builtin('subsasgn',w,S,val);
    case '{}'
        %w = tospace(w); % convert the internal representation of coefficients into space
        if length(S(1).subs)==1 && length(S)==1
            J = sort(w.depth);
            Jm = circshift(J(:),-1);Jm(1) = 0;
            M = 1 + sum(double(J(:)-Jm(:)).*(2.^(numel(J):-1:1)'-1));clear Jm;
            k = S(1).subs{1};
            if k<1||k>M
                error('WAVELETClass:subsref',...
                    'Not a supported subscripted reference');
            end
            if k<1||k>M
                error('WAVELETClass:subsref',...
                    'Not a supported subscripted reference');
            end
            if k==M
                w.coarse{end} = val;return;
            else
                counter = 1;
                for i = 1:length(w.wav)
                    for m = 1:length(w.wav{i})
                        if counter==k
                            w.wav{i}{m} = val;return;
                        end
                        counter = counter+1;
                    end
                end
            end
        elseif S(1).subs{1} == max(w.depth)+1,
            if ~all(size(val)==size(w.coarse{end}))
                error('WAVELETClass:subsref',...
                    'wrong dimensions');
            end
            w.coarse{end} = val;return;
        elseif length(S(1).subs)>1
            if ~all(size(val)==size(w.wav{S(1).subs{1}}{S(1).subs{2}}))
                error('WAVELETClass:subsref',...
                    'wrong dimensions');
            end
            w.wav{S(1).subs{1}}{S(1).subs{2}} = val;return;
        elseif length(S)>1
            if ~all(size(val)==size(w.wav{S(1).subs{1}}{S(2).subs{1}}))
                error('WAVELETClass:subsref',...
                    'wrong dimensions');
            end
            w.wav{S(1).subs{1}}{S(2).subs{1}} = val;return;
        else
            error('WAVELETClass:subsref',...
                'Not a supported subscripted reference');
        end
    otherwise
        error('WAVELETClass:subsref',...
            'Not a supported subscripted reference');
end