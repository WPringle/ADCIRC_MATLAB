function [Dr,Ds] = Dmatrices2D(N, r, s, V)
% Taken from JSH

[Vr, Vs] = GradVandermonde2D( N, r, s ) ;
Dr = Vr/V ; % Dr V = Vr
Ds = Vs/V ; % Ds V = Vs

return ;