% Eluer forward time stepping for surface eq
% ================ Variables ==================
% H_old -- Thickness at H^n 
% intv -- or vbar, average x velocity in vertical direction, equiv to qx/H
% dvdx -- d(vbar)/dx, upwinding based on the sign of qx. 
%               If qx>0, then dintvdx=(q(i)-q(i-1))/dx; if qx<0, then dintvdx=(q(i+1)-q(i))/dx;
% as -- accumulation on the top surface
% dt -- time step 
% ============== Return Values ================
% H_new -- new solution of thickness at H^{n+1}




function H_new = EulerImplicitThick(H_old, dx, intv, as, dt)
    % System size
    N = max(size(intv));
    % Assemble Matrix and rhs
    [Stiff, Force] = AssembleMatrix(dx, intv, as);
    
    % Solve for H_t+ Stiff*H =Force 
    Q = speye(N)/dt + Stiff;
    rhs = 1/dt*H_old + Force;
    H_new = Q \ rhs;
end