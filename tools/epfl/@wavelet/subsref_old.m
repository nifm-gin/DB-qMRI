function out = subsref(w, S)

depth = numel(S);

% detect the call of 'w(:)'
ToVector = (depth==1&&strcmp(S(1).type,'()')&&numel(S(1).subs)==1&&(strcmp(S(1).subs(1),':')));

if ToVector % provide the wavelet coefficients as a 1D vector
    w = tospace(w); % convert the internal representation of coefficients into space
    coeff = zeros(prod(w.size),1);
    Nold = 0;
    tmp = w.coarse{end}(:);N = Nold+length(tmp);
    coeff(Nold+1:N) = tmp; Nold = N;
    for i = 1:length(w.wav)
        for m = 1:length(w.wav{i})
            tmp = w.wav{i}{m}(:);N = Nold+length(w.wav{i}{m}(:));
            coeff(Nold+1:N) = tmp; Nold = N;
        end
    end
    return;
elseif strcmp(S(1).type,'{}')
    S(1).type
    S(1).subs
    if isa(S(1).subs{1},'double')
        if S(1).subs{1} == w.depth+1,out = w.coarse;else out = w.wav;end
        S = [S(1), S];
    elseif isa(S(1).subs,'char')
        if strcmp(S(1).subs,'end'), disp('coucou');out = w.coarse;end
    else
        error('bad call');
    end
    depth = depth+1;
elseif ~strcmp(S(1).type,'.')||~ischar(S(1).subs)
    error('bad call');
else
    switch S(1).subs
        case 'wav'
            out = w.wav;
        case 'coarse'
            out = w.coarse;
        case 'depth'
            out = w.depth;
            if depth==2&&strcmp(S(2).type,'()')
                out = subsref(out, S(2));
            end
            return;
        case 'size'
            out = w.size;
            if depth==2&&strcmp(S(2).type,'()')
                out = subsref(out, S(2));
            end
            return;
        case 'fourier'
            out = w.fourier;
            return;
        otherwise
            S(1).subs;
    end
end
w = tospace(w); % convert the internal representation of coefficients into space
for i = 2:depth
    S(i).type
    switch S(i).type
        case '()'
            tmp = ifftn(out);
            out = subsref(tmp, S(end)); % access pixels of the wavelet coeffs in the spatial domain
        case '{}'
            S(i).subs
            if numel(S(i).subs)~=1
                error('bad call');
            end
            if (i==2) % concerning the scale level
                out = out{S(i).subs{1}};
                if isa(out,'double') % case of coarse coeff
                    return;
                end
            elseif (i==3) % concerning the wavelet subband
                out = out{S(i).subs{1}}; % access wavelet subband in the spatial domain
            else
                error('bad call');
            end
        otherwise
            error('bad call');
    end
end

end