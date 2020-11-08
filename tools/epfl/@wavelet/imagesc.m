function imagesc(w)

%w = tospace(w);

switch numel(w.size)
    case 1
        error('1D plot to be defined!');
    case 2
        wavplot = zeros([w.size, 3]);
        indx = w.size(1);
        indy = w.size(2);
        MJ = max(w.depth);
        mJ = min(w.depth);
        for j = 1:MJ
            num = numel(w.wav{j});
            for s = 1:num
                tmp = w.wav{j}{s};
                if j<=mJ
                    switch s
                        case 1, x = 1:indx/2;y = indy/2+1:indy;
                        case 2, x = indx/2+1:indx;y = 1:indy/2;
                        case 3, x = indx/2+1:indx;y = indy/2+1:indy;
                    end
                else
                    error('still to program: case where decomposition depths are different depending on dimensions');
                end
                for c = 1:3
                    %size(wavplot(x,y,c))
                    %size(tmp)
                    wavplot(x,y,c) = tmp/max(abs(tmp(:))+eps);%
                end
                m(j,s) = max(abs(tmp(:)));
            end
            if j<=w.depth(1)
                wavplot(indx/2+1,1:indy,1) = NaN; % for the delimiting lines
                wavplot(indx/2+1,1:indy,2:3) = 0;
            end
            if j<=w.depth(2)
                wavplot(1:indx,indy/2+1,1) = NaN; % for the delimiting lines
                wavplot(1:indx,indy/2+1,2:3) = 0;
            end
            if j<=w.depth(2), indy = indy/2;end
            if j<=w.depth(1), indx = indx/2;end
        end
        tmp = w.coarse{end};
        x = 1:indx;y = 1:indy;
        for c = 1:3
            wavplot(x,y,c) = tmp/max(abs(tmp(:))+eps);%
        end
        m = [m(:); max(abs(tmp(:)))];
        wavplot = abs(wavplot);%/(max(m)+eps);%
        wavplot(find(isnan(wavplot))) = 1;
        imagesc(wavplot);title(['wavelet coefficients (modulus) rescaled per subband. max value: ' num2str(max(m)) '']);axis off;axis image;
    otherwise
        error('multidimensional plot to be defined!');
end

