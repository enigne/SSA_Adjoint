% D0 Operator of first derivative for N-vectors with space size dx 
% Central differences qx = (q(i+1)-q(i-1))/2dx
% BCFlag = 0 then B.C. is Direchlet(Default), so D(end,end) === 0
% BCFlag = 1 then B.C. is Neumann 
% BCFlag = 2 then B.C. is Periodical, so D(end,1) === 1
%
function D = Dcd(N, dx, BCFlag)
    % Set default value for dx and BCFlag if not supply
    switch nargin
        case 1
            dx = 1;
            BCFlag = 0;
        case 2
            BCFlag = 0;
    end

    % Central differences
    I = 0.5*speye(N);
    D = circshift(I,[0,1])-circshift(I,[1,0]);
    D(1,end) = 0;
    D(end,1) = 0;
    
    % Boundary Condition
    if BCFlag == 0 % Direhelet 
        D(1,1) = -1;
        D(1,2) = 1;
        D(end,end) = 1;
        D(end,end-1) = -1;    
    elseif BCFlag == 1 % Neumann 
        D(1,1) = 1;
        D(end,end) = 1;
        D(1,2) = 0;
        D(end,end-1) = 0;    
    elseif BCFlag == 2 % Periodical x(1)=x(N)
        D(1,end) = -0.5;
        D(end,1) = 0.5;
    end
    D = sparse(D)/dx;
end

