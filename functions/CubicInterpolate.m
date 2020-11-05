function y = CubicInterpolate( v_0, v_1, v_2, v_3, x )

N_0 = v_1;
N_1 = v_2;
G_0 = 0.5*( v_2 - v_0 );
G_1 = 0.5*( v_3 - v_1 );

A = N_0;
B = G_0;
C = -3*N_0 + 3*N_1 -2*G_0 - G_1;
D = 2*N_0 - 2*N_1 + G_0 + G_1;

x_2 = x*x;
x_3 = x_2 * x;

y = A + B*x + C*x_2 + D*x_3;