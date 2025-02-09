function P = Update_P(X, R, P, param)
gamma = param.gamma;
tau = param.tau;
U = L21Parameter(P',1e-5);
P = (tau*R*X') / (tau*X*X'+2*gamma*U);
