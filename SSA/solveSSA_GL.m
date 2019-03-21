function [u, H, beta, bxc, A, rhog, dx, n, m] = solveSSA_GL(N, b, dt, T)
    %%
    Nx = N + 1;
    L = 18e5;
    dx = L/(Nx-1);
    x = linspace(-dx,L,Nx)';
    
    %%
    if nargin < 2
        b =  -(-720 +778.5*(x/7.5e5) );
        if nargin < 3
            dt = 1;
            if nargin < 4
                T = 10000;
            end
        end
    end

    %% symmetric ice divide (extra grid point)
    b(1)=b(3); 
    bxc = diff(b)/dx;
    
    H=zeros(Nx,1)+10;
    
    %% Constants
    secperyear=365.25*24*3600;

    %% Physical Constants
    n = 3.0;
    rhoi = 900.0/(secperyear^2);
    rhow = 1000.0/(secperyear^2);
    g = 9.81*secperyear^2;
    m = 1./3.; % basal sliding exponent
    C = 7.624e6; % basal friction coefficient
    as = 0.3;
    rhog = rhow * g;
    A = 4.6416e-25*2; % Pa s^-1

    %% Unit convertion
    C = (C/(secperyear^m));
    A = A*secperyear;

    %% Solve SSA GL problem
    [~, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C, m, n, rhoi, rhow, g, as, dt, T);
    
    u = u(2:end);
    % H on stagger grid
    H = (H(1:end-1)+H(2:end)) * 0.5;
    beta = beta(2:end);
end

