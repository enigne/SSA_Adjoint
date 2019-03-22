clear
close all
%%
addpath('SSA')
%% Load final H and u from init file
load('DATA/SSAinit_N400.mat')

%% Setup restart
N_restart = 100;
H_mat = zeros(length(H), N_restart);
u_mat = zeros(length(u), N_restart);
gpos_vec = zeros(1, N_restart);
dt = 1;

%% Solve SSA GL problem
for i = 1: N_restart
    [gpos, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C, m, n, rhoi, rhow, g, as, dt, dt, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
    gpos_vec(i) = gpos;
end
%% For adjSSA you need the input
rhoig = rhoi*g;
sigma= 0.5e4; 
ist = 150;
n=3;
x = x(2:end);
Nx = length(x)+1;
Ascal = 1e-6;
% beta = beta(2:end);

%% Solve backward in time
psi_old = zeros(Nx-1, 1);
I = eye(Nx-1);

for i =  N_restart:-1:1
    % get the forward solutions
    u = u_mat(:,1);
    H = H_mat(:,1);
    u = u(2:end);
    % H on stagger grid
    H = (H(1:end-1)+H(2:end)) * 0.5;
    % construct Adjoint matrices
    [A11, A12, A21, A22, F1, F2]=constrauctAdjSSA(Nx-1,n,ist,sigma,u,H,mean(bxc),A,rhoig,dx);
    % Time stepping
    Q = [A11 - 1./dt .* I,	A12;
        A21,                A22;];
    rhs = [F1 - 1./dt .* psi_old;
        F2;];
    
    psifi = Q\rhs;
    psi_old = psifi(1:Nx-1);
    fi_old = psifi(Nx:2*Nx-2);    
    
    
    
    figure(1)
    subplot(2,1,1)
    plot(x, psi_old)
    subplot(2,1,2)
    plot(x, fi_old)
end

%%
% [psi,fi,wght,bwght,A,f,fic,eta,beta,fix,psix,hx,ux,B,G,K]=adjssa(Nx-1,n,ist,sigma,u,H,mean(bxc),A,rhoig,dx);
