function W = Update_W(R, W, A)
rho = 1;
rate_rho = 1.2;
iter = 1;
error = 1;
previousW = W;
Jw = W;
Zw = zeros(size(Jw));
while(error>1e-8&&iter<100)
    W = (R*A' + rho*Jw - rho*Zw) / (A*A' + rho*eye(size(A*A')));
    Jw = normcol_lessequal(W + Zw);
    Zw = Zw + W + Jw;
    rho = rate_rho*rho;
    error = mean(mean((previousW - W).^2));
    previousW = W;
    iter = iter + 1;
end