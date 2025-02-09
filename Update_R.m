function R = Update_R(X, Ys, Yt, R, P ,T, D, W, A, H, param)
alpha = param.alpha;
beta = param.beta;
lambda = param.lambda;
tau = param.tau;
DictSize = param.DictSize;
%% Obtain M
% Estimate sigma
C = length(unique(Ys));
Ns = size(Ys,2);
Nt = size(Yt,2);
Rs = R(:,1:Ns)';
Rt = R(:,(Ns+1):end)';
list_adist_c = [];
epsilon = 1e-3;
% A distance between marginal distribution
adist_m = adist(Rs, Rt);
% A distance between conditional distribution
for i=1:C
    index_i = Ys == i; 
    Rsi = Rs(index_i,:); 
    index_j = Yt == i; 
    Rtj = Rt(index_j,:); 
    adist_i = adist(Rsi,Rtj); 
    list_adist_c = [list_adist_c;adist_i];
end
adist_c = mean(list_adist_c);
sigma = adist_m / (adist_m + adist_c);
if sigma > 1
    sigma = 1;
elseif sigma <= epsilon
    sigma = 0;
end

% Compute M
% Compute marginal distribution M0
e = [1 / Ns * ones(Ns,1); -1 / Nt * ones(Nt,1)];
M0 = e * e' * C;
% Compute conditional distribution Mc
Mc = 0;
for c = reshape(unique(Ys),1,C)
    e = zeros(Ns + Nt,1);
    e(Ys == c) = 1 / length(find(Ys == c));
    e(Ns + find(Yt == c)) = -1 / length(find(Yt == c));
    e(isinf(e)) = 0;
    Mc = Mc + e * e';
end
M = sigma*M0 + (1-sigma)*Mc;

%% Update R
DictLabel = [];
for j=1:C
    DL(1,1:DictSize) = j;
    DictLabel = [DictLabel, DL];
end
L = Construct_L(D', DictLabel);
Temp_R = R;
Temp_RN = D'*T*X+tau*P*X+alpha*W*A+lambda*H.*R;
Temp_RD = D'*D*R+(tau+alpha+lambda)*R+beta*R*M+lambda*L*R;
R = Temp_R .* (Temp_RN./Temp_RD);




