var C K I N Y G;
predetermined_variables K;
parameters alpha beta gamma varphi rho delta;
varexo g_shock; 

alpha = 1/3;
beta = 0.99;
gamma = 1;
varphi = 1;
rho = 0.0;
delta = 0.02;
model;

C^(-gamma) = beta*C(+1)^(-gamma)*(1-delta+alpha*(K(+1)/N(+1))^(alpha-1));
C^(-gamma)*(1-alpha)*(K/N)^(alpha) = N^varphi;
//N = 1;
I + G + C = Y;
I = K(+1) - (1-delta)*K;
G = rho*G(-1) + g_shock;
Y = K^alpha*N^(1-alpha);

end;

initval;
C = 1;
N = 1;
K = 1;
G = 0;
I = 0.4;
Y = 1.2;
end;

steady(solve_algo = 3);

shocks;
var g_shock;
stderr  1;
end;
check;
stoch_simul(irf=10,order=1) C I Y G N; 
