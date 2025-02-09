function[R, P, Q, T, D, W, V, LM] = Initial_ADADPL(X, H, A, C, param)
DictSize = param.DictSize;
ReducedDim = param.ReducedDim;
gamma = param.gamma;

% Initialize synthetical dictionary D
rng(1,'twister');
D = normcol_equal(randn(ReducedDim, DictSize*C));

% Initialize analytical dictionary P
rng(2,'twister');
P = normcol_equal(randn(size(X,1), DictSize*C))';

% Initialize coding coefficients R
R = P*X;

% Initialize variable matrix Q
Q = H;

% Initialize embedding matrix W
rng(3,'twister');
W = normcol_equal(randn(size(R,1),size(A,1)));

% Initialize transformation matrix T
LM = zeros(ReducedDim, ReducedDim); % Lagrange multiplier
V = eye(size(X'*X))-1/size(X,2)*ones(size(X'*X));
temp_A = LM;
temp_B = (X*X'+gamma*eye(size(X*X')))/(X*V*X');
temp_C = (D*R*X')/(X*V*X');
T = sylvester(temp_A,temp_B,temp_C);
LM = LM + T*X*V*X'*T' - eye(size(LM));

