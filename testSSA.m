clear
close all
%%
addpath('SSA')
%% Load final H and u from init file
load('DATA/SSAinit_N1600.mat')
%% Setup restart
N_restart = 1;
uObs = 0;
HObs = 1 - uObs;
transientFlag = 0;
dt = 1;
%% Initialization
H_mat = zeros(length(H), N_restart);
u_mat = zeros(length(u), N_restart);
gpos_vec = zeros(1, N_restart);
%% Hack for m=1
C = C.*abs(u).^(m-1);
m = 1;
%% Solve SSA GL problem
for i = 1: N_restart
    [glInd, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C, m, n, rhoi, rhow, g, as, dt, dt, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
    gpos_vec(i) = glInd;
end
%% For adjSSA you need the input
rhoig = rhoi*g;
sigma= 1e2; 
ist = 700;
n = 3;

%% Solve backward in time
psi_old = zeros(N+1, 1);
phi_old = zeros(N+1, 1);
% Artificial viscosity
epsilon = 1e3;
Ascale = 1e0;
for i =  N_restart:-1:1
    % get the forward solutions
    u = u_mat(:,i);
    H = H_mat(:,i);
    glInd = gpos_vec(i);
    % set for adjoint
    xAdj = x(2:glInd);
    Nx = length(xAdj)+1;
    I = eye(Nx-1);
    %
    u = u(2:glInd);
    % H on stagger grid
    H = (H(1:glInd-1)+H(2:glInd)) * 0.5;
    % Beta
    Cbeta = beta(2:glInd);
    
    % construct Adjoint matrices
    [A11, A12, A21, A22, F1, F2, ux, eta, bb]=constrauctAdjSSAMatrices(...
        Nx-1,n,ist,sigma,u,H,mean(bxc),A,rhoig,dx, uObs, HObs, glInd, epsilon, Cbeta, m, 0, transientFlag);
    % Time stepping
    Q = [A11*Ascale - transientFlag*Ascale./dt .* I,	A12;
        A21*Ascale,                A22;];
    rhs = [F1 - transientFlag*Ascale./dt .* psi_old(1:glInd-1);
        F2;];

    psifi = Q\rhs;
    psi_old = psifi(1:Nx-1)*Ascale;
    phi_old = psifi(Nx:2*Nx-2);    
   
    wght = -phi_old .* u.^m;
    bwght = (Dm(Nx-1, dx)*psi_old).*u + (Dm(Nx-1,dx)* phi_old) .*eta .* ux+ rhoig*phi_old.*(Dm(Nx-1,dx)*H + bxc(1:Nx-1));
%     bwght = rhoig*(H.*(Dm(Nx-1,dx)*phi_old) + phi_old.*(Dm(Nx-1,dx)*H + bxc(1:Nx-1))) +F1;
        
    figure
    subplot(4,1,1)
    plot(xAdj, psi_old)
    ylabel('$\psi$', 'Interpreter','latex')
    subplot(4,1,2)
    plot(xAdj, phi_old)
    ylabel('$\phi$', 'Interpreter','latex')
    subplot(4,1,3)
    plot(xAdj, wght)
    ylabel('$\delta\beta$', 'Interpreter','latex')
    subplot(4,1,4)
    plot(xAdj, bwght )
end

%%
pert = sum(wght.*C(1:glInd-1)*0.01)*dx;

%%

t1 = rhoig*(H.*(Dm(Nx-1,dx)*phi_old));
t2 = rhoig* phi_old.*(Dm(Nx-1,dx)*H + bxc(1:Nx-1)) ;
t3 = F1;
% t4 = rhoig*phi_old.*(bxc(1:Nx-1));
figure
% range1 = [-0.8e-7,0.8e-7];
% range1 = [-0.2e-5,0.2e-5];
subplot(3,1,1)
plot(t1)
% ylim(range1)
subplot(3,1,2)
plot(t2)
% ylim(range1)
subplot(3,1,3)
plot(t3)
% ylim(range1)
% subplot(4,1,4)
% plot(t4)
% ylim(range1)
%% b weights
% 
% t1 = (Dp(Nx-1, dx)* psi_old) .*u;
% t2 = (Dp(Nx-1, dx)* phi_old) .*eta .* ux;
% t3 = rhoig*phi_old.*(Dp(Nx-1,dx)*H + bxc(1:Nx-1));
% figure
% subplot(3,1,1)
% plot(t1)
% subplot(3,1,2)
% plot(t2)
% subplot(3,1,3)
% plot(t3)
%%
% range2 = [-0.5e-4,0.5e-4];
% range2 = [-2e-6,2e-6];
% 
% t1 = Dp(Nx-1,dx)*(1./n.*H.*eta.* (Dm(Nx-1,dx)*phi_old));
% t2 =-m*bb.*phi_old;
% t3 = -H.*(Dp(Nx-1,dx)*psi_old);
% figure
% subplot(3,1,1)
% plot(t1)
% ylim(range2)
% subplot(3,1,2)
% plot(t2)
% ylim(range2)
% subplot(3,1,3)
% plot(t3)
% ylim(range2)
% 
% %%
% figure
% % plot(u./eta./H.*(Dp(Nx-1,dx)*(eta.*H) ))
% hold on
% % plot(eta.*ux)
% % plot(-rhog.*H)
% % hold on
% plot(1./(-rhog.*H))
% plot(phi_old)
