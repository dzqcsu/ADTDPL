function [Ts, Tt, D_Mat, Ps_Mat, Pt_Mat, TO] = TrainADTDPL2(Xs, Xt, Ys, y_t, Attribute, param)
X = blkdiag(Xs, Xt);
C = length(unique(Ys)); % number of class
%% Initialize
[R, P, Q, D, W, V, LM, H, A, Yt, T] = Initial_ADADPL(Xs, Xt, Ys, y_t, Attribute, param);
trustable = 0;
ds = size(Xs,1);
DictSize = param.DictSize;
TO = 0;
%% Update
for i=1:param.Iter
    tic
    % Update R
    [R, M, L] = Update_R(X, Ys, Yt, R, P ,T, D, W, A, H, param);
    % Obtain indicators
    [ind1,ind2] = Obtain_Ind(X, P, M, H, Q, L);
    % Update pseudo labels
    if sum(trustable)<length(y_t)
        [Yt, H, A, R, trustable] = Update_Label(R, Ys, y_t, H, A, Attribute, param, i, ind1, ind2);
    end
    % Update P
    P = Update_P(X, R, P, param);
    % Update Q
    Q = Update_Q(R, Q, H);
    % Update T
    [T, LM] = Update_T(X, R, D, V, LM, param);
    % Update D
    D = Update_D(X, R, T, D);
    % Update W
    W = Update_W(R, W, A);
    TT = toc;
    TO = TO+TT;
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
