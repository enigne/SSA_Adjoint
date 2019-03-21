function [u,res,beta] = SFMSolve(u,H,grl,taud,txx,A,C,gridx,N,m,n,D1m, subBeta, tol, matIter)
switch nargin
    case 13
        tol = 1e-6;
        maxIter = 5;
    case 14
        tol = min(tol,1e-6);
        maxIter = 5;
    case 15
        maxIter = max(5, maxIter);
end

u_old = u;
%%
res = zeros(maxIter,1);
for i=1:maxIter
    [u,beta] = SFMVelocity(u_old, H, -taud, N, A, C, m, n, gridx, grl, txx(end), D1m, subBeta);
    res(i) = norm(u_old - u)/norm(u_old);
    if (res(i) <= tol)
        break;
    end
    u_old = u;
end
