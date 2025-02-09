function[Ts, Tt, D_Mat, Ps_Mat, Pt_Mat] = TrainADTDPL1(Xs, Xt, Ys, Yt, Attribute, param)
%% Input
X = blkdiag(Xs, Xt);
C = length(unique(Ys)); % number of class
ds = size(Xs,1);
DictSize = param.DictSize;
Iter = param.Iter;
Hs = [];
Ht = [];
As = [];
At = [];
for i=1:C
    inds = find(Ys==i);
    indt = find(Yt==i);
    ns = length(inds);
    nt = length(indt);
    temps = ones(DictSize, ns);
    tempt = ones(DictSize, nt);
    Hs = blkdiag(Hs, temps);
    Ht = blkdiag(Ht, tempt);
    As(:,inds) = repmat(Attribute(i,:)',1,ns);
    At(:,indt) = repmat(Attribute(i,:)',1,nt);
end
H = [Hs, Ht];
A = [As, At];

%% Initialize
[R, P, Q, T, D, W, V, LM] = Initial_ADADPL(X, H, A, C, param);

%% Update
for i=1:Iter
    % Update R
    R = Update_R(X, Ys, Yt, R, P ,T, D, W, A, H, param);
    % Update P
    P = Update_P(X, R, P, param);
    % Update Q
    Q = Update_Q(R, Q, H);
    % Update T
    [T, LM] = Update_T(X, R, D, T, V, LM, param);
    % Update D
    D = Update_D(X, R, T, D);
    % Update W
    W = Update_W(R, W, A);
end

%% Output
Ts = T(:,1:ds);
Tt = T(:,(ds+1):end);
for j = 1:C
    ind = [((j-1)*DictSize+1):(j*DictSize)];
    D_Mat{j} = D(:,ind);
    Temp_P = P(ind,:);
    Ps_Mat{j} = Temp_P(:,1:ds);
    Pt_Mat{j} = Temp_P(:,(ds+1):end);
end
