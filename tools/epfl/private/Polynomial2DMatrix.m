%% Polynomial2DMatrix.m
%
% Function that provides the matrixt hat links the vector of parameters s
% to the values of D-degree polynomial at the points defined by R.
%
% INPUT:    * R: (2,m) matrix that defines m 2D points of evaluation
%           * D: degree of polynomial
%
% OUTPUT:   * M: (m,N) matrix that links the vector of parameters s that
%               defines a polynomial of degree D, to its values a the n
%               points.
%
% Note: In 2D we have N = (D+1)*(D+2)/2.
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 31-10-2009 (dd-mm-yyyy)

function M = Polynomial2DMatrix(R,D)

N = (D+1)*(D+2)/2;      % Number of parameters
m = size(R,2);          % Number of points

% construct the polynomial degrees for each parameter
or = zeros(2,N);
for d = 0:D
    n = d*(d+1)/2;
    for i = 0:d
        or(:,n+i+1) = [i;d-i];
    end
end

% build matrix
M = zeros(m,N);
for n = 1:N
    M(:,n) = R(1,:).^or(1,n).*R(2,:).^or(2,n);%prod(R(:,j).^or(:,n));
end