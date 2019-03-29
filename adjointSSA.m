% Adjoint SSA Solver with surface h coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'N'                 - Number of nodes;
%   'transientFlag'     - 1: trasient adjoint, 0: steady state;
%   'uObs'              - Flag to take observation on u;
%   'HObs'              - Flag to take observation on H;
%   'sigma'          	- width of the observation function F, it is 
%                           automaticaly scaled to such that the area is 1;
%   'epsilon'          	- Artificial viscosity for the right boundary close
%                           to grounding line;
% The return values:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2019-03-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function adjointSSA(N, transientFlag, uObs, HObs, sigma, epsilon)
if nargin < 6
    % Artificial viscosity
    epsilon = 1e3;
    if nargin < 5
        % width of the observation
        sigma= 0.5e4;
        if nargin < 4
            if nargin < 3
                uObs = 1;
                if nargin < 2
                    transientFlag = 0;
                end
            end
            HObs = 1 - uObs;
        end
    end
end
%%
addpath('SSA')
%% Load final H and u from init file
load(['DATA/SSAinit_N', num2str(N), '.mat'])
saveFlag = 1;
%% Setup restart
ist = [1:GLpos-1];
N_restart = 1;
dt = 1;
%% Initialization
Nist = length(ist);
H_mat = zeros(length(H), N_restart);
u_mat = zeros(length(u), N_restart);
gpos_vec = zeros(N_restart,1);
psi_mat = zeros(N+1, Nist, N_restart);
phi_mat = zeros(N+1, Nist, N_restart);
wght_mat = zeros(N+1, Nist, N_restart);
bwght_mat = zeros(N+1, Nist, N_restart);

%% Solve SSA GL problem
for i = 1: N_restart
    [glInd, H, u, ~]=FlowlineSSA(H, b, x, dx, Nx, A, C, m, n, rhoi, rhow, g, as, dt, dt, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
    gpos_vec(i) = glInd;
end
%% For adjSSA you need the input
rhoig = rhoi*g;
n = 3;
%% Solve backward in time
psi_old = zeros(N+1, 1);

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
        % forward solution
        u = u(2:glInd);
        % H on stagger grid
        H = (H(1:glInd-1)+H(2:glInd)) * 0.5;
        % construct Adjoint matrices
        [A11, A12, A21, A22, F1, F2, ux, eta]=constrauctAdjSSAMatrices(...
            Nx-1, n, ist(i), sigma, u, H, mean(bxc), A, rhoig, dx, uObs, HObs, glInd, epsilon);
        % Time stepping
        Q = [A11 - transientFlag*1./dt .* I,	A12;
            A21,                                A22;];
        rhs = [F1 - transientFlag*1./dt .* psi_old(1:glInd-1);
            F2;];
        % solve
        psifi = Q\rhs;
        psi = psifi(1:Nx-1);
        phi = psifi(Nx:2*Nx-2);
        
        % sensitivity
        wght = -phi .* u.^m;
        bwght = (Dp(Nx-1, dx)*psi).*u + (Dcd(Nx-1,dx)* phi) .*eta .* ux+ rhoig*phi.*(Dcd(Nx-1,dx)*H + bxc(1:Nx-1));
        
        psi_mat(1:glInd-1, i, t) = psi;
        phi_mat(1:glInd-1, i, t) = phi;
        wght_mat(1:glInd-1, i, t) = wght;
        bwght_mat(1:glInd-1, i, t) = bwght;
    end
end

%% save to data files
if saveFlag
    if uObs
        obsName = 'u';
    else
        obsName = 'H';
    end
    
    if transientFlag
        save(['DATA/SSAAdjoint_N', num2str(N), '_T', num2str(N_restart), '_', obsName ,'.mat'], ...
            'x', 'ist', 'psi_mat', 'phi_mat', 'wght_mat', 'bwght_mat');
    else
        save(['DATA/SSAAdjoint_N', num2str(N), '_',obsName , '.mat'], ...
            'x', 'ist', 'psi_mat', 'phi_mat', 'wght_mat', 'bwght_mat');
    end
end
