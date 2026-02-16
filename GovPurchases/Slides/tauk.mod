var C K I N Y tau_k;
predetermined_variables K;
parameters alpha beta gamma varphi rho delta disut;
varexo tau_k_shock; 

alpha = 0.7;
beta = 0.999;
gamma = 100000;
varphi = 1;
rho = 0.0;
delta = 0.01;
disut = 1;
model;

C^(-gamma) = beta*C(+1)^(-gamma)*(1-delta+(1-tau_k(+1))*alpha*(K(+1)/N(+1))^(alpha-1));
C^(-gamma) = disut*N^varphi;
//N = 1;
I + C = Y;
I = K(+1) - (1-delta)*K;
tau_k = tau_k_shock + tau_k(-1);
Y = K^alpha*N^(1-alpha);

end;

initval;
C = 1;
N = 1;
K = 6;
tau_k = 0;
I = delta*K;
Y = K^alpha*N^(1-alpha);
end;

steady(solve_algo = 3);

shocks;
var tau_k_shock;
stderr  0.01;
end;
check;
stoch_simul(irf=1000,order=1); 
