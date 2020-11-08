function bool = eq(A,B)
    size = prod(double(A.size~=B.size));
    depth = prod(double(A.J~=B.J));
    fam = ~strcmp(A.family,B.family);
    decim = xor(A.decimation,B.decimation);
    direc = xor(A.adjoint,A.inverse)~=xor(B.adjoint,B.inverse);
    if size||depth||fam||decim||direc
        bool = false;
    else
        bool = true;
    end
end