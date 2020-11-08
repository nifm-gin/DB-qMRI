% make.m

currentdir = pwd;
if isempty(strfind(currentdir,'private'))
    cd private/;
end
%CCFLAGS="-march=core2 -O2 -msse -msse2"
mex nuft_forw.cpp
mex nuft_back.cpp
mex nuft_kern.cpp
mex nuft_forw_gridding.cpp
mex nuft_back_gridding.cpp
mex myerfzparts_c.cpp
mex insidepoly_dblengine.c
cd(currentdir);

cd waveletstuff;

mex DownSampleDim.cpp
mex FilterSep.cpp
mex UpSampleDim.cpp
mex wavadd.cpp
mex wavmult.cpp
mex wavtovect.cpp

cd(currentdir);
