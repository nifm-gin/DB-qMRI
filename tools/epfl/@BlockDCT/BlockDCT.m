% Class of object that represents block discrete cosine transform.
% The classical operators for matrix multiplication are overloaded to
% permit some manipulations as for matrices.
%
% The class relies on the Matlab implementations.
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 10-08-2009 (dd-mm-yyyy)

classdef BlockDCT
    properties
        inverse = false;
        blocksize = 8;
        D;
    end
    
    methods
        function  A = BlockDCT(blksize)    % constructor
            A.inverse = false;
            A.blocksize = blksize;
            A.D = dctmtx(A.blocksize);
        end
        function s = display(A)
            if A.inverse == 0 %A*x
                direction = 'direct';
            else %At*x
                direction = 'inverse';
            end
            s = [direction ' ' num2str(A.blocksize) 'x' num2str(A.blocksize) ' block DCT'];
            disp(s);
        end
    end % methods
end % classdef
