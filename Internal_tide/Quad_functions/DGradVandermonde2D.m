function [V2Drr,V2Dss,V2Drs] = DGradVandermonde2D(N, r, s)
%
% Second derivative Vandermonde matrices 
%
%  Modal: \psi = 'Modal basis'
%
%     V2Drr = d^{2} ( \psi_{j} )/dr^{2}|_{ \boldsymbol{r}_{i} }
%     V2Dss = d^{2} ( \psi_{j} )/ds^{2}|_{ \boldsymbol{r}_{i} }
%     V2Drs = d^{2} ( \psi_{j} )/drds|_{ \boldsymbol{r}_{i} }
% DW
%
V2Drr = zeros(length(r),(N+1)*(N+2)/2) ;
V2Dss = zeros(length(r),(N+1)*(N+2)/2) ;
V2Drs = zeros(length(r),(N+1)*(N+2)/2) ; 

sk = 1 ;
for i = 0: N 
    for j = 0: N - i
        for l = 1: length(r)
            %
            [V2Drr(l,sk),V2Dss(l,sk),V2Drs(l,sk)] = DGradSimplex2DP(r(l), s(l), i, j) ;
            %
        end
        sk = sk + 1 ;
    end
end

return ;