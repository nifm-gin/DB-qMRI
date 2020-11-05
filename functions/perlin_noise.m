% Perlin noise algorithm
function I = perlin_noise( nbRandomPts, scale )

% First, generate random numbers
randoms = rand( nbRandomPts, 1 );

step  = 1/scale;
nbPts = scale*nbRandomPts;
I = zeros( nbPts, 1 );

for i = 1 : nbRandomPts

    v_1 = randoms( i );

    if( i == 1 )
        v_0 = randoms( nbRandomPts );
        v_2 = randoms( i+1 );
        v_3 = randoms( i+2 );
    elseif( i == nbRandomPts )
        v_0 = randoms( i-1 );
        v_2 = randoms( 1 );
        v_3 = randoms( 2 );
    elseif( i == nbRandomPts - 1 )
        v_0 = randoms( i-1 );
        v_2 = randoms( i+1 );
        v_3 = randoms( 1 );
    else
        v_0 = randoms( i-1 );
        v_2 = randoms( i+1 );
        v_3 = randoms( i+2 );
    end

    for j = 1 : scale
        I( scale*(i-1) + j, 1 ) = CubicInterpolate( v_0, v_1, v_2, v_3, j*step );
    end
end 