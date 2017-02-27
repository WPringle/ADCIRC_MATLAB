function [r,s] = xytors( x, y )
% Taken form JSH
% Given point in the equidistant triangle ==> 
% v1 = (-1,-1/sqrt(3)), v2 = (1,-1/sqrt(3)), v3 = (0,2/sqrt(3))
% corresponding points in the reference triangle
%  - use Barycnetric trnasformation
%      [1        2][L_{1}] = [x + 1]
%      [sqrt(3)  0][L_{3}]   [y + 1/sqrt(3)] 
% 
L1 = (sqrt(3.0)*y + 1.0)/3 ;
L2 = (-3.0*x - sqrt(3.0)*y + 2.0)/6.0 ; 
L3 = ( 3.0*x - sqrt(3.0)*y + 2.0)/6.0 ; 

r = -L2 + L3 - L1 ;
s = -L2 - L3 + L1 ; 

return ;