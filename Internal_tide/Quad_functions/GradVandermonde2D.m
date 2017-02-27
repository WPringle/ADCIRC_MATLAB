function [V2Dr,V2Ds] = GradVandermonde2D(N, r, s)
%
% Taken from JSH 
%
%  Modal: \psi = 'Modal basis'
%
%     V2Dr = d ( \psi_{j} )/dr|_{ \boldsymbol{r}_{i} }
%     V2Ds = d ( \psi_{j} )/ds|_{ \boldsymbol{r}_{i}}
%
V2Dr = zeros(length(r),(N+1)*(N+2)/2) ;
V2Ds = zeros(length(r),(N+1)*(N+2)/2) ;

sk = 1 ;
for i = 0: N 
    for j = 0: N - i
        for l = 1: length(r)
             [V2Dr(l,sk),V2Ds(l,sk)] = GradSimplex2DP(r(l), s(l), i, j) ;
        end
        sk = sk + 1 ;
    end
end

return ;