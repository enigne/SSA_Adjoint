clear
close all
%%
addpath('SSA')
%% Load final H and u from init file
load('DATA/SSAinit_N800.mat')

%% Setup restart
N_restart = 10;
H_mat = zeros(length(H), N_restart);
u_mat = zeros(length(u), N_restart);
gpos_vec = zeros(1, N_restart);
dt = 0.1;

%% Solve SSA GL problem
for i = 1: N_restart
    [glInd, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C, m, n, rhoi, rhow, g, as, dt, dt, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
    gpos_vec(i) = glInd;
end
%% For adjSSA you need the input
rhoig = rhoi*g;
sigma= 0.5e2; 
ist = 150;
n=3;
x = x(2:end);
Nx = length(x)+1;
Ascal = 1e-6;
% beta = beta(2:end);

%% Solve backward in time
psi_old = zeros(Nx-1, 1);
I = eye(Nx-1);
% Artificial viscosity
epsilon = 0;

for i =  N_restart:-1:1
    % get the forward solutions
    u = u_mat(:,i);
    H = H_mat(:,i);
    u = u(2:end);
    % H on stagger grid
    H = (H(1:end-1)+H(2:end)) * 0.5;
    % construct Adjoint matrices
%     [A11, A12, A21, A22, F1, F2]=constrauctAdjSSA(Nx-1,n,ist,sigma,u,H,mean(bxc),A,rhoig,dx);
    [A11, A12, A21, A22, F1, F2, ux, eta, beta]=constrauctAdjSSAMatrices(Nx-1,n,ist,sigma,u,H,mean(bxc),A,rhoig,dx, glInd, epsilon);
    % Time stepping
    Q = [A11 - 1./dt .* I,	A12;
        A21,                A22;];
    rhs = [F1 - 1./dt .* psi_old;
        F2;];
    
    psifi = Q\rhs;
    psi_old = psifi(1:Nx-1);
    fi_old = psifi(Nx:2*Nx-2);    
    
    psi_old(end) = 0;
    
    wght = -fi_old .* u.^m;
    
    figure(1)
    subplot(3,1,1)
    plot(x, psi_old)
    ylabel('$\psi$', 'Interpreter','latex')
    subplot(3,1,2)
    plot(x, fi_old)
    ylabel('$\phi$', 'Interpreter','latex')
    subplot(3,1,3)
    plot(x, wght)
    ylabel('$\delta\beta$', 'Interpreter','latex')
end

%%
pert = sum(wght.*C(1:end-1)*0.01)*dx;

%%
figure
t1 = Dcd(Nx-1,dx)*(1./n.*H.*eta.* (Dcd(Nx-1,dx)*fi_old));
t2 =-m*beta.*fi_old;
t3 = -H.*(Dcd(Nx-1,dx)*psi_old);
% subplot(3,1,1)
plot(t1+t2+t3)

% hold on
% plot(F2)
% plot(A21*psi_old+A22*fi_old -F2)
