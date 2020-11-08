function y = mtimes(A,x)

if 0%(length(A.htld_hat)==2)
    if ~xor(A.inverse == 0,A.adjoint == 0)
        y.coarse{1} = x;
        for i = 1:max(A.J)
            [y.coarse{i+1},y.wav{i}{1},y.wav{i}{2},y.wav{i}{3}] = dwt2(y.coarse{i},conj(A.htld_hat).',conj(A.gtld_hat).','mode','per');
        end
        y = wavelet(y,false);
    else % DWT inv
        for i = max(A.J):-1:1
            x.coarse{i} = idwt2(x.coarse{i+1},x.wav{i}{1},x.wav{i}{2},x.wav{i}{3},A.h_hat,A.g_hat,'mode','per');
        end
        y = x.coarse{1};
    end
else
    if A.inverse == 0 % DWT
        if A.adjoint == 0 %A*x
            x = fftn(x); % assumes the given object is represented in spatial domain
            y = tospace(wavelet(wav_anal(x,A.htldstr_hat,A.gtldstr_hat,A.J,A.decimation),true));
        else %At*x
            x = tofourier(x); % the wavelet coefficients should be represented in Fourier for the transform
            y = ifftn(wav_rec(x,A.htld_hat,A.gtld_hat,A.J,A.decimation));
        end
    else % DWT inv
        if A.adjoint == 0 %A^-1*x
            x = tofourier(x); % the wavelet coefficients should be represented in Fourier for the transform
            y = ifftn(wav_rec(x,A.h_hat,A.g_hat,A.J,A.decimation));
        else %A^-t*x
            x = fftn(x); % assumes the given object is represented in spatial domain
            y = tospace(wavelet(wav_anal(x,A.hstr_hat,A.gstr_hat,A.J,A.decimation),true));
        end
    end
end

function x_wav_hat = wav_anal(x_hat,low,high,J,decimation)

% Decomposition depth
jmax = size(low, 1);

x_wav_hat.coarse{1} = x_hat;

for j = 1:jmax
    [x_wav_hat.coarse{j+1}, x_wav_hat.wav{j}] = WavAnalysisSep(x_wav_hat.coarse{j}, low(j, :), high(j, :), decimation);
end

function x_hat = wav_rec(x_wav_hat,low,high,J,decimation)

% Decomposition depth
jmax = size(low, 1);

for j = jmax:-1:1
    x_wav_hat.coarse{j} = WavSynthesisSep(x_wav_hat.coarse{j+1}, x_wav_hat.wav{j}, low(j, :), high(j, :), decimation);
end
x_hat = x_wav_hat.coarse{1};