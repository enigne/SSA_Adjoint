% Assemble Matrices and Forcing term




function [Stiff, Force] = AssembleMatrix(dx, intv, as)
    % System size
    N = max(size(intv));
   
    % Determine the upwinding direction	
    upFlag = intv < 0;    % dhdx > 0;
    D1up = Dup(N, 1, upFlag, 1); 
    
    D = D1up/dx;    
    Stiff = spvardiag(D*intv);
    Force = as;
end