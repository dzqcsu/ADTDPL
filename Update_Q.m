function Q = Update_Q(R, Q, H)
G = 2*H - ones(size(H));
Q = 2*R - G.*Q;