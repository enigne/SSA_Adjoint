% Create spares diagonal matrix with vector v on the diagonal
% i > 0 upper triangle
% i < 0 lower triangle
% N size of targeting matrix
function Diag = spvaroffdiag(v,i,N)
    v = v(1:N-abs(i));
    if i > 0 
        rows = 1:N-i;
        columns = (i+1):N;
        Diag = sparse(rows,columns,v,N,N);
    elseif i < 0
        rows = (1-i):N;
        columns = 1:N+i;
        Diag = sparse(rows,columns,v,N,N);
    else 
        Diag = spvardiag(v);
    end
end