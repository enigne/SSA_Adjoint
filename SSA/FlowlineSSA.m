function [gpos, H, u, beta]=FlowlineSSA(H, b, x, dx, N, A, C, m, n, rhoi, rhow, g, as, dt, T)

% D operators
D1p = Dp(N);
D1m = Dm(N);


% Initialize 
u=zeros(N,1);
hstag=zeros(N,1);
xstag=zeros(N,1);
bstag=zeros(N,1);

% Constants
sea_level=0; % height of eustatic sea level

% time stepping
time_lapse = round(T/dt)+1;

for time_count = 1:time_lapse    
    % floating condition for ice sheet geometry determination
    haf=b-sea_level+H*rhoi/rhow; % height above floating
    hb=b;
    hb(haf<0)=sea_level-rhoi*H(haf<0)/rhow;
    s=hb+H;
    
    % geometry variables on u-grid
    for j=1:N-1
        %         hstag(j,1)=(H(j+1)+H(j))/2.; % ice thickness on u-grid
        %         xstag(j,1)=(x(j)+x(j))/2.; % distance on u-grid
        %         bstag(j,1)=(b(j+1)+b(j))/2.;
        hstag(j,1)=(H(j+1)+H(j))/2.; % ice thickness on u-grid
        xstag(j,1)=(x(j+1)+x(j))/2.; % distance on u-grid
        bstag(j,1)=(b(j+1)+b(j))/2.;
    end
    hstag(N)=NaN;
    xstag(N)=NaN;
    bstag(N)=NaN;
    
    % floating condition on u-grid
    haf=bstag-sea_level+hstag*rhoi/rhow; % height above floating
    grl=ones(N,1)*2; % initialize grl-vector
    grl(haf>0)=0; % grounded
    grl(haf<0)=2; % floating
    
    % Subgrid determination of grounding line position on u-grid
    grlj=N;
    for j=1:N-2
        if grl(j)<.5 && grl(j+1)>1.5
            grl(j)=1; % last grounded grid point grlj
            grlj=j;
        end
    end
    
    % Longitudinal and driving stresses
    slope=diff(s)/dx; % slope on u-grid
    slope(N)=NaN;
    taud = -rhoi*g.*hstag.*slope; % driving stress on u-grid
    txx = 0.25*rhoi*g*H*(1.-rhoi/rhow); % on h-grid
    
    if time_count==1 % initialization of velocity field with SIA
        ub=C.^(-1./m)*abs(taud.^(1./m-1.)).*taud; % basal velocity (m/year) on u-grid
        u=ub+2./(n+2.)*A*hstag.*abs(taud).^(n-1.).*taud; % ice velocity in ice sheet
    end
    
    % grounding line index grlj
    f =(sea_level-bstag).*rhow./(rhoi*hstag);
    df = (f(grlj+1) - f(grlj))/dx;
    x_grl = (1-f(grlj)+df*x(grlj))./df;
    x_sub = x_grl - x(grlj);
    subBeta = dx/x_sub;
    if subBeta > 2.5
        subBeta = 0;
    end
    % Solve SSA
    [u,~,beta] = SFMSolve(u,H,grl,taud,txx,A,C,dx,N,m,n,D1m,subBeta);
    
    % Solve thickness Equation
    H = StaggerThick(H, dx, u, as, dt,grl, D1m, D1p);
    
end %% end of time stepping

disp(x(grlj)/1e3);
gpos=x(grlj);

% 
figure(1)
plot(x(2:end),H(2:end)+hb(2:end)); hold on;
plot(x(2:end),b(2:end),'linewidth',2);
plot(x(2:end),hb(2:end));
hold off;
ylim([-1000,5000])
xlabel('x')
ylabel('H')
figure
plot(x(2:end),u(2:end),'linewidth',2);
xlabel('x')
ylabel('v')
end
