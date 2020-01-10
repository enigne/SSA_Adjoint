% Adjoint SSA Solver with surface h coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'N'                 - Number of nodes;
%   'transientFlag'     - 1: trasient adjoint, 0: steady state;
%   'uObs'              - Flag to take observation on u;
%   'HObs'              - Flag to take observation on H;
%   'MacayealFalg'      - To use MacAyeal's formulation of the adjoint SSA
%   'ist'               - The index of the observation, can be a scalar or
%                           an array;
%   'T_final'          	- The final time of the simulation, the adjoint
%                           equation starts from that point;
%   'dtAdj'             - The time step;
%   'obsT'              - The length of the observation;
%   'seasonType'        - The type of the seasonal variations on C, 0- no,
%                           1-sine, 2-cosin, 3-square
%   'sigma'          	- width of the observation function F, it is
%                           automaticaly scaled to such that the area is 1;
%   'epsilon'          	- Artificial viscosity for the right boundary close
%                           to grounding line;
% The return values:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: Old: 2019-03-27
%       Updated: 2020-01-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function adjointSSA(N, transientFlag, uObs, HObs, MacayealFalg, ist, ...
    T_final, dtAdj, obsT, seasonType, amplitude, sigma, epsilon)
if nargin < 13
    % Artificial viscosity
    epsilon = 1e3;
    if nargin < 12
        % width of the observation
        sigma= 1e2;
        if nargin < 11
            amplitude = 0.5;
            if nargin < 10
                seasonType = 0;
                if nargin < 9
                    obsT = 0.1;
                    if nargin < 8
                        dtAdj = 1;
                        if nargin < 7
                            T_final = 1;
                            if nargin < 6
                                ist = 0;
                                if nargin < 5
                                    MacayealFalg = 0;
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
                        end
                    end
                end
            end
        end
    end
end
%%
addpath('SSA')
%% Load final H and u from init file
load(['DATA/SSAinit_N', num2str(N), '.mat'])
saveFlag = 1;
%% Pertubation in time
N_restart = ceil(T_final / dtAdj);
sInd = [0: (N_restart-1)];

if seasonType == 0
    season = 0 * sInd;
elseif seasonType == 1
    % sine [-1,1]
    season = sin(sInd*dtAdj*2*pi);
elseif seasonType == 2
    % cosine [-1,1]
    season = cos(sInd*dtAdj*2*pi);
elseif seasonType == 3
    % square [0, 1]
    season = (sin(sInd*dtAdj*2*pi) > 0)*2-1;
end

%% Hack for m=1
C = C.*abs(u).^(m-1);
dC = amplitude.*C;
m = 1;

%% For adjSSA you need the input
rhoig = rhoi*g;
n = 3;

%% Initialization
Nist = length(ist);

H_mat = zeros(length(H), N_restart);
u_mat = zeros(length(u), N_restart);
beta_mat = zeros(length(C), N_restart);
gpos_vec = zeros(N_restart,1);

psi_mat = zeros(N+1, N_restart, Nist);
phi_mat = zeros(N+1, N_restart, Nist);
wght_mat = zeros(N+1, N_restart, Nist);
bwght_mat = zeros(N+1, N_restart, Nist);

if MacayealFalg
    v_mat = zeros(N+1, N_restart, Nist);
    vWght_mat = zeros(N+1, N_restart, Nist);
end


%% Solve SSA GL problem
for i = 1: N_restart
    [glInd, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C+dC*season(i), 1, n, rhoi, ...
        rhow, g, as, dtAdj, dtAdj, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
    beta_mat(:, i) = beta;
    gpos_vec(i) = glInd;
end


%% Solve backward in time
psi_old = zeros(N+1, 1);

for i =  1:Nist
    for t =  N_restart:-1:1
        % get the forward solutions
        u = u_mat(:,t);
        H = H_mat(:,t);
        glInd = gpos_vec(t);
        Cbeta = beta_mat(2:glInd, t);
        
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
            Nx-1, n, ist(i), sigma, u, H, mean(bxc), A, rhoig, dx, uObs, ...
            HObs, glInd, epsilon, Cbeta, m, MacayealFalg, transientFlag);
        
        % Time dependent observation
        if t <= (N_restart - obsT/dtAdj)
            F1 = 0*F1;
            F2 = 0*F2;
        end
        
        % Time stepping
        Q = [A11 - transientFlag*1./dtAdj .* I,	A12;
            A21,                                A22;];
        rhs = [F1 - transientFlag*1./dtAdj .* psi_old(1:glInd-1);
            F2;];
        
        % solve
        psifi = Q\rhs;
        psi = psifi(1:Nx-1);
        phi = psifi(Nx:2*Nx-2);
        
        % sensitivity
        wght = -phi .* u.^m;
        bwght = (Dm(Nx-1, dx)*psi).*u + (Dm(Nx-1,dx)* phi) .*eta .* ux + ...
            rhoig*phi.*(Dm(Nx-1,dx)*H + bxc(1:Nx-1));
        
        % save sensitivities
        psi_mat(1:glInd-1, t, i) = psi;
        psi_old = psi;
        phi_mat(1:glInd-1, t, i) = phi;
        wght_mat(1:glInd-1, t, i) = wght;
        bwght_mat(1:glInd-1, t, i) = bwght;
        
        % solve steady state SSA using Macayeal's formulation
        % only valid for u observation
        if MacayealFalg
            v = A22 \ F2;
            v_mat(1:glInd-1, t, i) = v;
            wght = -v .* u.^m;
            vWght_mat(1:glInd-1, t, i) = wght;
        end
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
        if MacayealFalg
            save(['DATA/SSA_Macayeal_Adjoint_N', num2str(N), '_T', num2str(N_restart),...
                '_', obsName, '_S', num2str(seasonType), '.mat'], 'x', 'ist', 'dtAdj', ...
                'N_restart', 'psi_mat', 'phi_mat', 'wght_mat', 'bwght_mat', ...
                'v_mat', 'vWght_mat', 'beta_mat');
            
        else
            save(['DATA/SSAAdjoint_N', num2str(N), '_T', num2str(N_restart), '_', obsName , ...
                '_S', num2str(seasonType), '_A', num2str(floor(amplitude*100)), '.mat'], ...
                'x', 'ist', 'dtAdj', 'N_restart', 'psi_mat', 'phi_mat', ...
                'wght_mat', 'bwght_mat', 'beta_mat', 'season');
        end
        
    else
        if MacayealFalg
            save(['DATA/SSA_Macayeal_Adjoint_N', num2str(N), '_', obsName , ...
                '.mat'], 'x', 'ist', 'dt', 'N_restart', 'psi_mat', 'phi_mat', ...
                'wght_mat', 'bwght_mat', 'v_mat', 'vWght_mat', 'beta_mat');
        else
            save(['DATA/SSAAdjoint_N', num2str(N), '_',obsName , '.mat'], ...
                'x', 'ist', 'dt', 'N_restart', 'psi_mat', 'phi_mat', ...
                'wght_mat', 'bwght_mat', 'beta_mat' );
        end
    end
end
