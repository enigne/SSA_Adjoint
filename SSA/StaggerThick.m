% Stagger grid for thickness eq
% ================ Variables ==================
% H_old -- Thickness at H^n 
% intv -- or vbar, average x velocity in vertical direction, equiv to qx/H
% dvdx -- d(vbar)/dx, upwinding based on the sign of qx. 
%               If qx>0, then dintvdx=(q(i)-q(i-1))/dx; if qx<0, then dintvdx=(q(i+1)-q(i))/dx;
% as -- accumulation on the top surface
% dt -- time step 
% ============== Return Values ================
% H_new -- new solution of thickness at H^{n+1}




function H_new = StaggerThick(H_old, dx, intv, as, dt, grl, D1m, D1p)
    % System size
    N = max(size(intv));
    % Assemble Matrix and rhs
%     D1p = Dp(N);
%     D1m = Dm(N);
    intvm = circshift(intv,[1,1]); % no matter intv is row or column vector

    uDiag = spvardiag(intv);    
    uDiagm = spvardiag(intvm);
    
    % Solve for H_t+ Stiff*H =Force 
    Q = speye(N) + (D1m*uDiag+D1p*uDiagm)*dt/2/dx;
    rhs = H_old + as*dt;
    
    % Boundary conditions at H(1)
    Q(1,1) = 1;
    Q(1,2) = 0; 
    Q(1,3) = -1;
    rhs(1) = 0;
    
    % Boundary conditions at H(N)
    if grl(N)>1.5
        %H(N)=H(N-1); % ice shelf
        Q(N,N) = 1;
        Q(N,N-1) = -1;
        rhs(N) = 0;
    else
        %H(N)=0; % ice sheet
        Q(N,N) = 1;
        Q(N,N-1) = 0;
        rhs(N) = 0;
    end
    
    P = 1+intv-intvm;
    P(1) = 1;
    P(N) = 1;
    P = spvardiag(P);
    Q = P*Q;
    rhs = P*rhs;
    
    H_new = Q \ rhs;
%     H_new(1) = H_new(3);
end