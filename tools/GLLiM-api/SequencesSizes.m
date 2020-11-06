function [Sequences,s1,s2,t,slices] = SequencesSizes(Sequences)

% Fabien Boux - 11/2020

narginchk(1, 1);

% Identify the MRI data sizes
switch length(size(Sequences))
    case 4
        [s1,s2,t,slices] = size(Sequences);
    case 3
        [s1,s2,t]   = size(Sequences);
        slices      = 1;
    case 2
        [s1,t]      = size(Sequences);
        s2          = 1;
        Sequences   = reshape(Sequences, s1,s2,t);
        slices      = 1;
    otherwise
        error('Invalid Sequences argument size')
        
end