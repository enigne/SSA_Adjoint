% D+ Operator of first derivative for N-vectors with space size dx 
% one side finite differences qx = (q(i+1)-q(i))/dx
% BCFlag = 0 then B.C. is Direchlet(Default), so D(end,end) === 0
% BCFlag = 1 then B.C. is Neumann 
% BCFlag = 2 then B.C. is Periodical, so D(end,1) === 1
%
function D = Dp(N, dx, BCFlag)
    % Set default value for periodFlag
    switch nargin
        case 1
            dx = 1;
            BCFlag = 0;
        case 2
            BCFlag = 0;
    end
    
    % One side finite differences
    I = speye(N);
    D = -I + circshift(I,[0,1]);
    D(end,1) = 0;

    if BCFlag == 0      % Direhelet 
        D(end,end) = 1;
    elseif BCFlag == 1	% Neumann 
        D(end,end) = 1;
    elseif BCFlag == 2	% Periodical x(1)=x(N)
        D(end,1) = 1;
    end

	D = sparse(D);
    if dx ~= 1
        D = D./dx;
    end
end