% Create spares diagonal matrix with vector v on the diagonal
% 
function Diag = spvardiag(v)
    N = length(v);
    p = 1:N;
	Diag = sparse(p,p,v);
end