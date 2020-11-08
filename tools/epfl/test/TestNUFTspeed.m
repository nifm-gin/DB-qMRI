%% TestNUFTspeed.m
%
% Testing the speed of several implementations of the matrix-vector
% multiplications m = E*x, y = EH*m, and z = M*x, where E is a particular
% non uniform discrete Fourier transform matrix, ET is its Hermitian transpose
% and M = EH*E (with a specific implementation).
%
% Copyright, Matthieu Guerquin-Kern, 2012


%clear all;close all;clc;
disp('Performing the tests for non uniform Fourier transforms speed...');
FOV = .24;
methods = {'matlab script','DFT mex','nuft gg','nuft gg fast'};
repeat = 5;
dim = round(2.^(3:.5:6));
tE0 = zeros(numel(dim),numel(methods));
tEH0 = zeros(numel(dim),numel(methods));
tEHE0 = zeros(numel(dim),numel(methods));

for d = 1:numel(dim)
    mxsize = dim(d);
    x = phantom(mxsize);
    mx = size(x);
    
    w = GenerateSpiralTraj(FOV,FOV/mxsize,1,1,10,1);
    [k,mx] = TrajInGridUnits(w,FOV,mx);
    
    %% Computations
    

%   testing E0

    for i = 1:numel(methods)
        param = [];
        param.k = k;
        param.method = methods{i};
        [m,param] = E(x, param);
        for r = 1:repeat
            t0 = clock();
            [m,param] = E0(x, param);
            tE0(d,i) = tE0(d,i) + etime(clock(),t0);
        end
    end
    tE0(d,i) = tE0(d,i)/repeat;
    
    %%
%  testing EH0
    
    for i = 1:numel(methods)
        param = [];
        param.k = k;
        param.method = methods{i};
        param.res = mx;
        [y,param] = EH0(m(:), param);
        for r = 1:repeat
            t0 = clock();
            [y,param] = EH0(m(:), param);
            tEH0(d,i) = tEH0(d,i) + etime(clock(),t0);
        end
    end
    tEH0(d,i) = tEH0(d,i)/repeat;
    
    %%
%    testing EHE0
    
    for i = 1:numel(methods)
        param = [];
        param.k = k;
        param.method = methods{i};
        param.res = mx;
        [z,param] = EHE0(x, param);
        for r = 1:repeat
            t0 = clock();
            [z,param] = EHE0(x, param);
            tEHE0(d,i) = tEHE0(d,i) + etime(clock(),t0);
        end
    end
    tEHE0(d,i) = tEHE0(d,i)/repeat;
end
%%
figure;loglog(dim,tE0,'-*');title('computation time for E0');legend(methods,'Location','NorthWest');set(gca,'XTick',dim);axis([dim(1) dim(end) min(tE0(1,:)) max(tE0(end,:))]);
figure;loglog(dim,tEH0,'-*');title('computation time for EH0');legend(methods,'Location','NorthWest');set(gca,'XTick',dim);axis([dim(1) dim(end) min(tEH0(1,:)) max(tEH0(end,:))]);
figure;loglog(dim,tEHE0,'-*');title('computation time for EHE0');legend(methods,'Location','NorthWest');set(gca,'XTick',dim);axis([dim(1) dim(end) min(tEHE0(1,:)) max(tEHE0(end,:))]);
