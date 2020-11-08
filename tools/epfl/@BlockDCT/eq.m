function bool = eq(A,B)
    size = A.blocksize==B.blocksize;
    direc = A.inverse==B.inverse;
    if size&&direc
        bool = true;
    else
        bool = false;
    end
end