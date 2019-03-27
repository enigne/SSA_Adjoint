clear
close all
%%
addpath('SSA')
saveFlag = 1;
%%
N = 1600;
Nx = N + 2;
L = 16e5;
dx = L/(N);
x = linspace(-dx,L,Nx)';

%%
b =  -(-720 +778.5*(x/7.5e5));
dt = 1;
T = 30000;

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
A = 4.6416e-25; % Pa s^-1

%% Unit convertion
C = (C/(secperyear^m))*ones(size(x));
A = A*secperyear;

%% Solve SSA GL problem
[GLpos, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C, m, n, rhoi, rhow, g, as, dt, T);

%% Save final H and u to mat file
if saveFlag 
    save(['DATA/SSAinit_N', num2str(N) ,'.mat'])
end
