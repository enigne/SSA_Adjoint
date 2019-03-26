clear
close all
%%
addpath('SSA')
%% Load final H and u from init file
load('DATA/SSAinit_N400.mat')
saveFlag = 1;
%% Setup restart
ist = [1:260];
N_restart = 10;
transientFlag = 1;
H_mat = zeros(length(H), N_restart);
u_mat = zeros(length(u), N_restart);
dt = 0.1;

Nist = length(ist);
psi_mat = zeros(N+1, Nist, N_restart);
phi_mat = zeros(N+1, Nist, N_restart);
wght_mat = zeros(N+1, Nist, N_restart);
bwght_mat = zeros(N+1, Nist, N_restart);

%% Solve SSA GL problem
for i = 1: N_restart
    [glInd, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C, m, n, rhoi, rhow, g, as, dt, dt, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
end
%% For adjSSA you need the input
rhoig = rhoi*g;
sigma= 0.5e4;
n=3;
x = x(2:end);
Nx = length(x)+1;
epsilon = 1000;

%% Solve steady states adjoint SSA
psi_old = zeros(Nx-1, 1);
I = eye(Nx-1);
% get the forward solutions
u = u(2:end);
% H on stagger grid
H = (H(1:end-1)+H(2:end)) * 0.5;
for i =  1:Nist
    for t =  N_restart:-1:1

    % construct Adjoint matrices
    [A11, A12, A21, A22, F1, F2, ux, eta]=constrauctAdjSSAMatrices(Nx-1,n,ist(i),sigma,u,H,mean(bxc),A,rhoig,dx,glInd,epsilon);
    % Time stepping
    Q = [A11 - transientFlag*1./dt .* I,	A12;
        A21,                                A22;];
    rhs = [F1 - transientFlag*1./dt .* psi_old;
        F2;];
    
    psifi = Q\rhs;
    
    psi = psifi(1:Nx-1);
    phi = psifi(Nx:2*Nx-2);
    wght = -phi .* u.^m;
    bwght = (Dp(Nx-1, dx)*psi).*u + (Dcd(Nx-1,dx)* phi) .*eta .* ux+ rhoig*phi.*(Dcd(Nx-1,dx)*H + bxc);
    
    psi_mat(:, i, t) = psi;
    phi_mat(:, i, t) = phi;
    wght_mat(:, i, t) = wght;
    bwght_mat(:, i, t) = bwght;
    end
end

% %%
% figure
% subplot(2,2,1)
% plot(psi_mat)
% ylabel('$\psi$', 'Interpreter','latex')
% subplot(2,2,2)
% plot(phi_mat)
% ylabel('$\phi$', 'Interpreter','latex')
% subplot(2,2,3)
% plot(wght_mat)
% ylabel('$\delta C$ sensitivity', 'Interpreter','latex')
% subplot(2,2,4)
% plot(bwght_mat)
% ylabel('$\delta b$ sensitivity', 'Interpreter','latex')

%%
if saveFlag 
if transientFlag
    save(['DATA/SSAAdjointTransient_T', num2str(N_restart) ,'_N', num2str(N) ,'.mat'], 'x', 'ist', 'psi_mat', 'phi_mat', 'wght_mat', 'bwght_mat');
else
    save(['DATA/SSAAdjoint_N', num2str(N) ,'.mat'], 'x', 'ist', 'psi_mat', 'phi_mat', 'wght_mat', 'bwght_mat');
end
end
