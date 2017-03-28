function Out = nearestinterpOTPS(T)
    J = find(T(3,:) == -99999 ); 
    if ~isempty(J)
        K = find(T(3,:) ~= -99999 );
        IDX = knnsearch([T(2,K)',T(1,K)'],[T(2,J)',T(1,J)']);
        T(3,J) = T(3,K(IDX));
        T(4,J) = T(4,K(IDX));     
    end
    % Final output in column vector form
    Out = T(3:4,:)';
end

