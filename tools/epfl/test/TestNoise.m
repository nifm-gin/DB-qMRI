%% TestNoise.m
%
% A random correlation matrix is made. A noise draw is correlated accordingly.
% The evaluation of the correlation matrix out of the noise data is checked.
%
% Copyright, Matthieu Guerquin-Kern, 2012

disp('Performing the test for noise synthesis and estimation...');
Nc = 4;
Ndata = 1e6;
% Generate correlation matrix
Psi=randn(Nc,Nc)+1j*randn(Nc,Nc);
Psi = (Psi'+Psi)/2;
[U, S, V] = svd (Psi);
Psi = V*abs(S)*V';
Psi = (Psi+Psi')/2;
% Correlate a white Gaussian noise draw
sPsi = chol(Psi);
fprintf('the following number should be as close to 0 as possible: %g\n',abs(trace(Psi-sPsi'*sPsi)/trace(Psi)));
n = randn(Ndata,Nc)+1j*randn(Ndata,Nc);n = n./std(n(:));
n=n*sPsi;
fprintf('the following number should be quite close to 0: %g\n',abs(trace(Psi-n'*n/Ndata)/trace(Psi)));
% Estimate covariance from the correlated noise draw
V = EstimateCovarianceMatrix(n);
fprintf('the following number should be as close to 0 as possible: %g\n',abs(trace(V*Psi-eye(Nc))/trace(Psi)));
