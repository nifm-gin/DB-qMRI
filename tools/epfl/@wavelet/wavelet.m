% Class that represents the wavelet coefficients obtained after a wavelet
% transform of type DWT. The classical linear operations are overloaded to
% permit the manipulation of these objects as vectors.
% Internaly, the wavelet coefficients may be stored in the Fourier domain.
%
% See also: DWT
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 17-07-2009 (dd-mm-yyyy)

classdef wavelet
    properties
        size; % size of the underlying object (numel(size) is the number of dimensions)
        depth; % depth of the decomposition (for each dimension)
        coarse; % cell array of the different coarse approximations
        wav; % cell array of the different wavelet levels
        fourier; % Are the coefficients stored as Fourier coeffs?
        decimated; % Are the coefficients decimated?
    end
    
    methods
        function  w = wavelet(x,IsFourier)
            w.coarse = x.coarse;
            w.wav = x.wav;
            if nargin == 2% constructor from structure
                s = size(x.coarse{1});
                ind = find(s>1);
                w.size = s(ind);
                w.decimated = [any(size(x.coarse{end})<size(x.coarse{1})), any(size(x.coarse{end})<size(x.coarse{end}))];
                if all(w.decimated)
                    s = size(x.coarse{end});
                    ind = find(s>1);
                    w.depth = uint8(log2(w.size./s(ind)));
                else
                    w.depth = length(w.wav)*ones(1,numel(w.size)); % assuming same decomposition depth in all dimensions
                end
                w.fourier = IsFourier;
            else % constructor by copy
                w.size = x.size;
                w.depth = x.depth;
                w.fourier = x.fourier;
                w.decimated = x.decimated;
            end
        end
        
        function display(w)
            size_strg = [];
            for i = 1:length(w.size)
                size_strg = [size_strg 'x' num2str(w.size(i))];
            end
            depth_strg = [];
            for i = 1:length(w.depth)
                depth_strg = [depth_strg 'x' num2str(w.depth(i))];
            end
            disp(['wavelet coefficients of a ' num2str(length(w.size)) '-D object of dimensions ' size_strg(2:end) ' with decomposition depths ' num2str(depth_strg(2:end))]);
        end
        
        function w = tofourier(w) % convert the internal representation of coefficients into Fourier
            if w.fourier == false % represented in space
                MJ = length(w.wav);
                for j = 1:MJ
                    num = numel(w.wav{j});
                    for s = 1:num
                        w.wav{j}{s} = fftn(w.wav{j}{s});
                    end
                end
%                 for j = 1:MJ+1
%                     w.coarse{j} = fftn(w.coarse{j});
%                 end
                w.coarse{MJ+1} = fftn(w.coarse{MJ+1});
                w.fourier =  true;
            end
        end
        
        function w = tospace(w)% convert the internal representation of coefficients into space
            if w.fourier == true % represented in Fourier
                MJ = length(w.wav);
                for j = 1:MJ
                    num = numel(w.wav{j});
                    for s = 1:num
                        w.wav{j}{s} = ifftn(w.wav{j}{s});
                    end
                end
%                 for j = 1:MJ+1
%                     w.coarse{j} = ifftn(w.coarse{j});
%                 end
                w.coarse{MJ+1} = ifftn(w.coarse{MJ+1});
                w.fourier =  false;
            end
        end
    end
end