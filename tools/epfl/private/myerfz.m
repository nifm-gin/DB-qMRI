function e = myerfz(R,IMPLEMENTATION)
% This function computes e = Erfz(R+1j*R) for R real
%
% Based upon the Matlab implementation of the error function for complex
% numbers by Marcel Leutenegger ( (c) January 2008, LGPL version 2.1).
% https://documents.epfl.ch/users/l/le/leuteneg/www/MATLABToolbox/ErrorFunc
% tion.html
if isempty(R)
    disp('found empty input, returning the input');
    e=R;
elseif isreal(R)
    e=nan(size(R));
    n=isfinite(R); % case 0 (nan input): return nan;
    if any(n==0), disp('found nan');end
    k=(n == uint8(1)); % case 1 (finite input): return erf(R) + parts(R)
    e(k)=erf(R(k)); % compute the first term (matlab's erf for reals)
    k = k&(abs(R)>eps); % locate non-zero inputs
    if any(k(:))
        if nargin<2, IMPLEMENTATION = 'fastest';end
        if ~strcmp(IMPLEMENTATION,'mfun')
        e(k)=e(k) + myerfzparts(R(k),IMPLEMENTATION);
        else
            e(k)= myerfzparts(R(k),IMPLEMENTATION);
        end
        %n=(R < uint8(0)); % locate negative inputs
        %R(n)=conj(R(n));
    end
else
    e=nan(size(R));
    n=isfinite(R); % case 0 (nan input): return nan;
    if any(n==0), disp('found nan');end
    k=(n == uint8(1)); % case 1 (finite input): return erf(R) + parts(R)
    e(k)=erf(real(R(k))); % compute the first term (matlab's erf for reals)
    k = k&(abs(R)>eps); % locate non-zero inputs
    if any(k(:))
        if nargin<2, IMPLEMENTATION = 'fastest';end
        if ~strcmp(IMPLEMENTATION,'mfun')
        %n=(R < uint8(0)); % locate negative inputs
        e(k)=e(k) + myerfzparts(R(k),IMPLEMENTATION);
        else
            e(k)= myerfzparts(R(k),IMPLEMENTATION);
        end
        %R(n)=conj(R(nend
    end
end

function e=myerfzparts(R,IMPLEMENTATION)
% This function computes e = Erf(R+1j*R) - Erf(R) for R nonzero real
%
% Based upon the Matlab implementation of the error function for complex
% numbers by Marcel Leutenegger ((c) January 2008, LGPL version 2.1).
% https://documents.epfl.ch/users/l/le/leuteneg/www/MATLABToolbox/ErrorFunc
% tion.html

%fprintf('statistics: min R: %f, max R: %f, mean R: %f, std R: %f\n',min(R(:)),max(R(:)),mean(R(:)),std(R(:)));
%figure;semilogy(abs(R(:)),'.');title('modulus of R');drawnow;

if nargin<2, IMPLEMENTATION = 'fastest';end
switch IMPLEMENTATION
    case {'fastest','best'}
        e=myerfzparts_c(R);
    case 'devel'
        e=myerfzparts_dev_c(R);
    case 'mfun'
        % INSANELY SLOW for large arrays
        % Becomes unstable for some R>6 (not all)
        e=mfun('erf',R);
    case 'matlab'
        if isreal(R)
            R2=R.*R;
            I= R;
            aI = abs(R);
        else
            I = imag(R);
            aI = abs(I);
            R = real(R);
            R2 = R.*R;
        end
        e2iRI=exp(complex(0,-2*R.*aI));
        E=(1 - e2iRI)./(2*pi*R);
        E(~R)=0;
        F=0;
        Hr=0;
        Hi=0;
        N=sqrt(1 - 4*log(eps/2));
        for n=1:ceil(N)
            H=n*n/4;
            H=exp(-H)./(H + R2);
            F=F + H;
            H=H.*exp(-n*aI);
            Hi=Hi + n/2*H;
            Hr=Hr + H;
        end
        e=exp(-R2).*(E + R.*F/pi - e2iRI.*complex(R.*Hr,Hi)/(2*pi));
        clear('E','F','H*');
        R3=R2 + log(2*pi);
        Gr=0;
        Gi=0;
        M=2*aI;
        n=max(1,floor(M - N));
        M=ceil(max(M + N - n));
        for m=0:M
            n1=n/2;
            n2=n1.*n1;
            G=exp(n.*aI - n2 - R3 - log(n2 + R2));
            Gi=Gi - n1.*G;
            Gr=Gr + G;
            n=n + 1;
        end
        e=e - e2iRI.*complex(R.*Gr,Gi);
        n=R == uint8(0);
        e(n)=e(n) + complex(0,aI(n)./pi);
        n = (I<0);
        e(n) = conj(e(n));
end