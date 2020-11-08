%% TestNUFTaccuracy.m
%
% Testing the accuracy of several implementations of the matrix-vector
% multiplications m = E*x, y = EH*m, and z = M*x, where E is a particular
% non uniform discrete Fourier transform matrix, EH is its Hermitian transpose
% and M = EH*E (with a specific implementation).
% The reference implementation does not use gridding.
% We check the consistency of the operations using the property
%	||E*x||² = <EH*E*x,x> = <M*x,x>   >= 0.
%
% Copyright, Matthieu Guerquin-Kern, 2012


close all;
%clear all;clc;
disp('Performing the tests for non uniform Fourier transforms accuracy...');
FOV = .24;
mx = [64, 123];
DefineBrain;	% check that the code handles nonsquare images and odd dimensions
x = RasterizePhantom(Brain,mx);

w = GenerateSpiralTraj(FOV,FOV./mx,1,1,10,1);
%k=w*FOV/2/pi;
[k,mx] = TrajInGridUnits(w,FOV,mx);
%figure(98);plot(k(:,1),k(:,2),'.-');axis square;title(sprintf('Spiral traj: %d x %d',mx(1),mx(2)));

methods = {'matlab script','DFT mex','nuft gg','nuft gg fast'};

%% testing E0

m = cell(1,numel(methods));
for i = 1:numel(methods)
    param = [];
    param.k = k;
    param.method = methods{i};
    m{i} = E0(x, param);
end
SERdB = log10(zeros(1,numel(methods)));
for i = 2:numel(methods)
    SERdB(i) = -20*log10(norm(m{i}(:)-m{1}(:))/norm(m{1}(:)));
end
str = ' SER (dB) |';
ttle = sprintf('E0        | %.7s | %.7s  | %.7s  | %.7s  |\n______________________________________________________',methods{1},methods{2},methods{3},methods{4});
disp(ttle);
disp([str, num2str(SERdB)]);

%% testing EH0
y = cell(1,numel(methods));
for i = 1:numel(methods)
    param = [];
    param.k = k;
    param.method = methods{i};
    param.res = mx;
    y{i} = EH0(m{1}(:), param);
end
SERdB = log10(zeros(1,numel(methods)));
consistency_imag = zeros(1,numel(methods));
consistency_real = zeros(1,numel(methods));
tmp = y{1}(:)'*x(:);
ref = norm(m{i}(:))^2/prod(mx);
consistency_imag(1) = -10*log10(abs(imag(tmp))/ref);
consistency_real(1) = -10*log10(abs(ref-tmp)/ref);
for i = 2:numel(methods)
    SERdB(i) = -20*log10(norm(y{i}(:)-y{1}(:))/norm(y{1}(:)));
    tmp = y{i}(:)'*x(:);
    consistency_imag(i) = -10*log10(abs(imag(tmp))/ref);
    consistency_real(i) = -10*log10(abs(ref-tmp)/ref);
end
ttle = sprintf('EH0       | %.7s | %.7s  | %.7s  | %.7s  |\n______________________________________________________',methods{1},methods{2},methods{3},methods{4});
str = ['SER (dB) |';' Const 1 |';' Const 2 |'];
disp(' ');
disp(ttle);
disp([str, num2str([SERdB;consistency_imag;consistency_real])]);
%% testing EHE0
z = cell(1,numel(methods));
for i = 1:numel(methods)
    param = [];
    param.k = k;
    param.method = methods{i};
    param.res = mx;
    [z{i},param] = EHE0(x, param);
end
SERdB = zeros(1,numel(methods));
consistency_imag = zeros(1,numel(methods));
consistency_real = zeros(1,numel(methods));
for i = 1:numel(methods)
    SERdB(i) = -20*log10(norm(z{i}(:)-y{1}(:))/norm(y{1}(:)));
    tmp = z{i}(:)'*x(:);
    consistency_imag(i) = -10*log10(abs(imag(tmp))/ref);
    consistency_real(i) = -10*log10(abs(ref-tmp)/ref);
end
ttle = sprintf('EHE0     | %.7s | %.7s  | %.7s  | %.7s  |\n______________________________________________________',methods{1},methods{2},methods{3},methods{4});
disp(' ');
disp(ttle);
disp([str, num2str([SERdB;consistency_imag;consistency_real])]);
