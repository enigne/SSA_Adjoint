clear
% close all
%%
addpath('SSA')
N = 1600;

%% Load final H and u from init file
load(['DATA/SSAinit_N', num2str(N), '.mat'])
saveFlag = 0;

%% Setup restart
uObs = 0;
HObs = 1 - uObs;
transientFlag = 1;
MacayealFalg = 0;

%% For adjSSA you need the input
rhoig = rhoi*g;
sigma= 1e3;
ist = [450];
n = 3;
dt = 0.05/4;
T_final = 10.5;
obsT = 0.1; % time period of the observation, from final time backward

%% Pertubation in time
N_restart = (T_final / dt);
sInd = [0: (N_restart-1)];
seasonType = 1;

if seasonType == 0
    season = 0 * sInd;
elseif seasonType == 1
    % sine [-1,1]
    season = sin(sInd*dt*2*pi);
elseif seasonType == 2
    % cosine [-1,1]
    season = cos(sInd*dt*2*pi);
elseif seasonType == 3
    % square [0, 1]
    season = (sin(sInd*dt*2*pi) > 0)*2-1;
end

%% Hack for m=1
pertubation = 0.5;
lWin = x(1);
rWin = x(end);
C = C.*abs(u).^(m-1);
dC = pertubation.*C.*((x>=lWin)&(x<=rWin));
m = 1;

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

u_ref = u;
H_ref = H;

%% Solve SSA GL problem
for i = 1: N_restart
    [glInd, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C+dC*season(i), 1, n, rhoi, ...
        rhow, g, as, dt, dt, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
    beta_mat(:, i) = beta;
    gpos_vec(i) = glInd;
end

%% Solve backward in time
psi_old = zeros(N+1, 1);
phi_old = zeros(N+1, 1);
% Artificial viscosity
epsilon = 1e3;

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
        if t <= (N_restart - obsT/dt)
            F1 = 0*F1;
            F2 = 0*F2;
        end
        
        % discrete system with implicit Euler
        Q = [A11 - 1./dt .* I,      A12;
            A21,                    A22;];
        
        rhs = [F1 - 1./dt .* psi_old(1:glInd-1);
            F2;];
        
        % solve
        psifi = Q\rhs;
        psi = psifi(1:Nx-1);
        phi = psifi(Nx:2*Nx-2);
        
        % sensitivity
        wght = -phi .* u.^m;
        bwght = (Dm(Nx-1, dx)*psi).*u + (Dm(Nx-1,dx)* phi) .*eta .* ux+ rhoig*phi.*(Dm(Nx-1,dx)*H + bxc(1:Nx-1));
        
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
                '_', obsName, '_S', num2str(seasonType), '.mat'], 'x', 'ist', 'dt', ...
                'N_restart', 'psi_mat', 'phi_mat', 'wght_mat', 'bwght_mat', ...
                'v_mat', 'vWght_mat', 'beta_mat');
            
        else
            save(['DATA/SSAAdjoint_N', num2str(N), '_T', num2str(N_restart), '_', obsName , ...
                '_S', num2str(seasonType), '.mat'], 'x', 'ist', 'dt', 'N_restart', 'psi_mat',...
                'phi_mat', 'wght_mat', 'bwght_mat', 'beta_mat');
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


%% Plot
% close all

tMesh = linspace(0, N_restart*dt, N_restart);
[Xm, Tm] = meshgrid(x(2:end), tMesh);

figure
colormap(jet)

subplot(2, 2, 1)
imagesc(x(2:end), tMesh, psi_mat')
xlim([0, x(GLpos)])
xlabel('x')
ylabel('t')
title('psi')
colorbar

subplot(2, 2, 2)
imagesc(x(2:end), tMesh, ((phi_mat')))
xlim([0, x(GLpos)])
xlabel('x')
ylabel('t')
title('v')
colorbar

subplot(2, 2, 3)
imagesc(x(2:end), tMesh, wght_mat')
xlim([0, x(GLpos)])
xlabel('x')
ylabel('t')
title('dC weights')
colorbar

subplot(2, 2, 4)
imagesc(x(2:end), tMesh, bwght_mat')
xlim([0, x(GLpos)])
xlabel('x')
ylabel('t')
title('db weights')
colorbar

%% slice Plot along x=x(ist)
figure

subplot(2, 2, 1)
plot(tMesh, psi_mat(ist+1,:))
hold on
plot(tMesh, psi_mat(ist,:))
plot(tMesh, psi_mat(ist+2,:))
xlabel('t')
ylabel('psi')
legend({'x(j)', 'x(j-1)', 'x(j+1)'}, 'location', 'northwest')

subplot(2, 2, 2)
plot(tMesh, phi_mat(ist+1,:))
hold on
plot(tMesh, phi_mat(ist,:))
plot(tMesh, phi_mat(ist+2,:))
xlabel('t')
ylabel('v')
legend({'x(j)', 'x(j-1)', 'x(j+1)'}, 'location', 'northwest')

subplot(2, 2, 3)
plot(tMesh, wght_mat(ist+1,:))
hold on
plot(tMesh, wght_mat(ist,:))
plot(tMesh, wght_mat(ist+2,:))
xlabel('t')
ylabel('dC weights')
legend({'x(j)', 'x(j-1)', 'x(j+1)'}, 'location', 'northwest')

subplot(2, 2, 4)
plot(tMesh, bwght_mat(ist+1,:))
hold on
plot(tMesh, bwght_mat(ist,:))
plot(tMesh, bwght_mat(ist+2,:))
xlabel('t')
ylabel('db weights')
legend({'x(j)', 'x(j-1)', 'x(j+1)'}, 'location', 'northwest')

