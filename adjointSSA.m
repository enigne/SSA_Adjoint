clear
close all
%%
addpath('SSA')
%% Load final H and u from init file
load('DATA/SSAinit_N400.mat')
saveFlag = 0;
%% Setup restart
ist = [1:266];
Nist = length(ist);
psi_mat = zeros(N+1, Nist);
phi_mat = zeros(N+1, Nist);
wght_mat = zeros(N+1, Nist);
bwght_mat = zeros(N+1, Nist);

%% For adjSSA you need the input
rhoig = rhoi*g;
sigma= 0.5e4;
n=3;
x = x(2:end);
Nx = length(x)+1;

%% Solve steady states adjoint SSA
I = eye(Nx-1);
% get the forward solutions
u = u(2:end);
% H on stagger grid
H = (H(1:end-1)+H(2:end)) * 0.5;
for i =  1:Nist
    % construct Adjoint matrices
    [A11, A12, A21, A22, F1, F2, ux, eta]=constrauctAdjSSAMatrices(Nx-1,n,ist(i),sigma,u,H,mean(bxc),A,rhoig,dx);
    % Time stepping
    Q = [A11, A12;
        A21, A22;];
    rhs = [F1;
        F2;];
    
    psifi = Q\rhs;
    
    psi = psifi(1:Nx-1);
    phi = psifi(Nx:2*Nx-2);
    wght = -phi .* u.^m;
    bwght = (Dp(Nx-1, dx)*psi).*u + (Dcd(Nx-1,dx)* phi) .*eta .* ux+ rhoig*phi.*(Dcd(Nx-1,dx)*H + bxc);
    
    psi_mat(:, i) = psi;
    phi_mat(:, i) = phi;
    wght_mat(:, i) = wght;
    bwght_mat(:, i) = bwght;
end

%%
figure
subplot(2,2,1)
plot(psi_mat)
ylabel('$\psi$', 'Interpreter','latex')
subplot(2,2,2)
plot(phi_mat)
ylabel('$\phi$', 'Interpreter','latex')
subplot(2,2,3)
plot(wght_mat)
ylabel('$\delta C$ sensitivity', 'Interpreter','latex')
subplot(2,2,4)
plot(bwght_mat)
ylabel('$\delta b$ sensitivity', 'Interpreter','latex')

%%
if saveFlag 
    save(['DATA/SSAAdjoint_N', num2str(N) ,'.mat'], 'x', 'ist', 'psi_mat', 'phi_mat', 'wght_mat', 'bwght_mat');
end
