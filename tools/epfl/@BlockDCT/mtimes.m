function y = mtimes(A,x)

if size(x)~=A.blocksize*floor(size(x)/A.blocksize);
    error(['Block DCT: incompatible data size, should be factor of ' num2str(A.blocksize) ' in each dimension']);
end

if A.inverse==1
    y =  blockproc(x,[A.blocksize A.blocksize],@(block_struct) A.D'*block_struct.data*A.D);
else
    y =  blockproc(x,[A.blocksize A.blocksize],@(block_struct) A.D*block_struct.data*A.D');
end