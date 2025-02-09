function D = Update_D(X, R, T, D)
rho = 1;
rate_rho = 1.2;
iter = 1;
error = 1;
previousD = D;
Jd = D;
Zd = zeros(size(Jd));
while(error>1e-8&&iter<100)
    D = (T*X*R' + rho*Jd - rho*Zd) / (R*R' + rho*eye(size(R*R')));
    Jd = normcol_lessequal(D + Zd);
    Zd = Zd + D + Jd;
    rho = rate_rho*rho;
    error = mean(mean((previousD - D).^2));
    previousD = D;
    iter = iter + 1;
end