classdef finite_difference_operator
    % finite_difference_operator immitates a matrix that performes a finite
    % difference calculation.
    %   The constuctor takes an array of directions (e.g. [1, 2, 3]). When
    %   multiplying the operator (mtimes) to an image, the finite
    %   difference is calculated along the dimensions 1, 2 and 3 and the
    %   result of each is concatted along dimension 1.
    %   E.g. in order to calculate a gradient, the adjoint of the operator
    %   can be use by calling A'. The sum over all dimension is taken.
    %   This operator uses periodic boundary conditions.
    %
    %   Example: A = finite_difference_operator([1 2]);
    %            x = rand(256);
    %            d = A*x;  % the result has the dimension [2 256 256]
    %            g = A'*d; % the result has again the dimension [256 256]
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (c) Jakob Asslaender, August 2016
    % New York University School of Medicine, Center for Biomedical Imaging
    % University Medical Center Freiburg, Medical Physics
    % jakob.asslaender@nyumc.org
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        direction = []; % Dimenensions along which the operation is perfomed
        adjoint = 0;    % Boolean indicating if the matrix is adjoint
    end
    
    methods
        function  A = finite_difference_operator(direction)
            % Constructor:
            %   A = finite_difference_operator(direction)
            %   Input: direction = vector of integer: [1 2]
            %          Indicates that the finite difference is calculated
            %          along the dimensions 1 and 2 and the result is
            %          summed up
            %   Output: Object
            
            A.direction = direction;
        end
        
        function A = ctranspose(A)
            % Sets the objection in the adjoint mode.
            % Called when using A'
            A.adjoint = ~A.adjoint;
        end
        
        function Q = mtimes(A,B)
            % Performs the actual calculation of the finite difference.
            %   Called when using A*x or A'*x
            %   Input: A is usually the finite_difference_operator
            %          B is the image if A is not adjoint.
            %          B is the finite diff. of an image if A is adjoint
            %   Output: Finite difference or image itself, if A is not or
            %           is adjoint, respectively.
            
            if isa(A, 'finite_difference_operator')
                if A.adjoint==1
                    if isvector(B)
                        Q = B - circshift(B,-1);
                    else
                        s = size(B);
                        Q = zeros(s(2:end));
                        for id = 1:length(A.direction)
                            Q = Q + squeeze(B(id,:,:,:,:)) - circshift(squeeze(B(id,:,:,:,:)),-1,A.direction(id));
                        end
                    end
                else
                    if isvector(B)
                        Q = B - circshift(B,1);
                    else
                        Q = [];
                        for id = 1:length(A.direction)
                            Q = cat(1, Q, reshape(B - circshift(B,1,A.direction(id)), [1 size(B)]));
                        end
                    end
                end
            else   % now B is the operator and A is the vector
                Q = mtimes(B',A')';
            end
        end
    end
end