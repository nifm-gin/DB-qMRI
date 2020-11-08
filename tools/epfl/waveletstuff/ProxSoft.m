function w = ProxSoft(T, w)

w = sign(w) .* min(T/2, abs(w));