%	Schoof' full model velocity solver
%	Solves for (2) in Schoof, 2007, Ice sheet grounding line dynamics:
%   steady states, stability and hysteresis
%   Author: CG
%   Date: 2016-05-03

% u=shelfystream(u,h,grl,taud,txx,A,c,gridx,secperyear,maxx,m,n);
% u = SFMVelocity(u, h, taud, maxx, A, c, m, n, gridx, grl, txx);
%
function [u,beta] = SFMVelocity(u, H, rhs, N, A, C, m, n, dx, grl, txx, D1m, subBeta)

    %must be D1m for upwinding on right boundary
    uxn = (D1m*u/dx).^2;
    uxn(1)=uxn(3);
    uxn(N)=uxn(N-1);
%     uxn(uxn>0)=uxn(uxn>0).^((1-n)/n*0.5);
    uxn(uxn<=0) = 1;
    uxn=uxn.^((1-n)/n*0.5);
    uxn(uxn<1) = 1;

    etaH = 0.5*H*A^(-1/n).*uxn; % mu on h-grid viscosity
    etaHp = circshift(etaH,-1);
    
    beta = C.*abs(u).^(m-1);
    beta((grl>=2)) = 0;
    beta(1)=beta(2);
    gl = find(grl==1);
    
    if (subBeta)
       beta(gl+1) = beta(gl)/subBeta;
    end
    
    dn = 2*(etaH+etaHp)./(dx^2);
    cen = -2*dn - beta;
    up = dn;
    
    delta = 2*(etaHp - etaH)./(dx^2);
    dn = dn - delta;
    up = up + delta;
    
    % B.C. u(1)=-u(2)
    up(1) = 1;
    cen(1) = 1;
    rhs(1) = 0;
    % B.C. u(N) = u(N-1)+txx
    cen(N) = 1;
    dn(N) = -1;
    rhs(end) = (A*txx^n)*dx;
    
    %
    B = [cen,[0;up(1:N-1)],[dn(2:N);0]];
    index = [0;1;-1];
    sysQ = spdiags(B,index,N,N);

    P = spvardiag(1.0./cen);
    sysQ = P * sysQ;
    rhs = P * rhs;
    u = sysQ\rhs;

end


