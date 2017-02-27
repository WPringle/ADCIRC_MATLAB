function [Drr,Dss,Drs] = D2matrices2D(N, r, s, V)
% 
%  2-nd differnation matrix
%
%  by DW
%  Note: V -- Vander 
[Vrr, Vss, Vrs] = DGradVandermonde2D( N, r, s ) ;

Drr = Vrr/V ; % /V % D2r V = Vr
Dss = Vss/V ; % /V % D2s V = Vs
Drs = Vrs/V ; % /V %

return ;