function w = SoftThresh(w,T)

if any(T<0)
    error('the threshold value should be positive!');
end

if isa(w,'wavelet')
    w = w - ProxSoftThresh(w,T);
else
    w = w - ProxSoft(T,w);
end