function  res = ctranspose(A)

A.inverse = xor(A.inverse,1);
res = A;