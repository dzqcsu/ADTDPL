function L=Construct_L(X,gnd)
options = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.NeighborMode = 'Supervised';
options.gnd = gnd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.WeightMode = 'HeatKernel';
options.t = 10^0;   
W = constructW(X,options);
Z=diag(aa);
L=Z-W;
