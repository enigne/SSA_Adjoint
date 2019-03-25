% D1 Operator of first derivative for N-vectors with space size dx 
% upwind scheme: if q> 0 qx = (q(i)-q(i-1))/dx; if q<0, qx = (q(i+1)-q(i))/dx;
% upFlag = 0 for q > 0; upFlag = 1 for q < 0;
% BCFlag = 0 then B.C. is Direchlet(Default), so D(end,end) === 0
% BCFlag = 1 then B.C. is Neumann 
% BCFlag = 2 then B.C. is Periodical, so D(end,1) === 1
%
function D = Dup(N, dx, upFlag, BCFlag)
    % Set default value for BCFlag
    switch nargin
        case 3
            BCFlag = 0;
    end
    
    % Get the operator only
    D1p = Dp(N, 1, BCFlag);
    D1m = Dm(N, 1, BCFlag);
        
    upSign = (upFlag > 0);
    % Sign matrix
    Pup = spvardiag(upSign);
    Pdown = spvardiag(ones(N,1)- upSign);
 
    % Upwind shceme
    D = Pup*D1p + Pdown*D1m;
    
    % Scalling
    D = D/dx;
end

