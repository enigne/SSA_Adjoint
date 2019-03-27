clear
close all
%%
addpath('SSA')
%% Load final H and u from init file
load('DATA/SSAinit_N400.mat')
saveFlag = 1;
%% Setup restart
ist = [1:260];
N_restart = 1;
transientFlag = 1;
H_mat = zeros(length(H), N_restart);
u_mat = zeros(length(u), N_restart);
gpos_vec = zeros(N_restart,1);
dt = 1;

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
    gpos_vec(i) = glInd;
end
%% For adjSSA you need the input
rhoig = rhoi*g;
sigma= 0.5e2;
n=3;

%% Solve backward in time
psi_old = zeros(N+1, 1);
% Artificial viscosity
epsilon = 1e4;

for i =  1:Nist
    for t =  N_restart:-1:1
        % get the forward solutions
        u = u_mat(:,t);
        H = H_mat(:,t);
        glInd = gpos_vec(t);
        % set for adjoint
        xAdj = x(2:glInd);
        Nx = length(xAdj)+1;
        I = eye(Nx-1);
        
        u = u(2:glInd);
        % H on stagger grid
        H = (H(1:glInd-1)+H(2:glInd)) * 0.5;
        % construct Adjoint matrices
        [A11, A12, A21, A22, F1, F2, ux, eta]=constrauctAdjSSAMatrices(Nx-1,n,ist(i),sigma,u,H,mean(bxc),A,rhoig,dx,glInd,epsilon);
        % Time stepping
        Q = [A11 - transientFlag*1./dt .* I,	A12;
            A21,                                A22;];
        rhs = [F1 - transientFlag*1./dt .* psi_old(1:glInd-1);
            F2;];
        
        psifi = Q\rhs;
        
        psi = psifi(1:Nx-1);
        phi = psifi(Nx:2*Nx-2);

        phi(ist(i)+5:end)=0;

        wght = -phi .* u.^m;
        bwght = (Dp(Nx-1, dx)*psi).*u + (Dcd(Nx-1,dx)* phi) .*eta .* ux+ rhoig*phi.*(Dcd(Nx-1,dx)*H + bxc(1:Nx-1));
        
        psi_mat(1:glInd-1, i, t) = psi;
        phi_mat(1:glInd-1, i, t) = phi;
        wght_mat(1:glInd-1, i, t) = wght;
        bwght_mat(1:glInd-1, i, t) = bwght;
    end
end

%%
if saveFlag
    if transientFlag
        save(['DATA/SSAAdjointTransient_T', num2str(N_restart) ,'_N', num2str(N) ,'.mat'], 'x', 'ist', 'psi_mat', 'phi_mat', 'wght_mat', 'bwght_mat');
    else
        save(['DATA/SSAAdjoint_N', num2str(N) ,'.mat'], 'x', 'ist', 'psi_mat', 'phi_mat', 'wght_mat', 'bwght_mat');
    end
end
