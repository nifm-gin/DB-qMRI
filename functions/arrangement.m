function y=arrangement(v,n)
    m=length(v);
    y=zeros(m^n,n);
    for k = 1:n
        y(:,k) = repmat(reshape(repmat(v,m^(n-k),1),m*m^(n-k),1),m^(k-1),1);
    end
end