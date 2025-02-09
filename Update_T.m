function [T, LM] = Update_T(X, R, D, T, V, LM, param)
gamma = param.gamma;
% Update T
temp_A = LM;
temp_B = (X*X'+gamma*eye(size(X*X')))/(X*V*X');
temp_C = (D*R*X')/(X*V*X');
T = sylvester(temp_A,temp_B,temp_C);
% Update Lagrange multiplier
LM = LM + T*X*V*X'*T' - eye(size(LM));