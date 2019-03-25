% D2 Operator of second derivative for N-vectors with space size dx 
% BCFlag = 0 then B.C. is Direchlet(Default), so D(end,end) === 0
% BCFlag = 1 then B.C. is Neumann 
% BCFlag = 2 then B.C. is Periodical, so D(end,1) === 1
%
function D = D2(N, dx, BCFlag)
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
    D = circshift(I,[1,0]) - 2*I + circshift(I,[0,1]);
    D(1,1) = 1;
    D(end,end) = 1;
  
    if BCFlag == 0      % Direhelet 
%         D(1,1) = 1;
    elseif BCFlag == 1	% Neumann 
%         D(1,1) = 1;
    elseif BCFlag == 2	% Periodical x(1)=x(N)
        D(1,end) = -1;
    end

	D = sparse(D);
    if dx ~= 1
        D = D./dx./dx;
    end
end