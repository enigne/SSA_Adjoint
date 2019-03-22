function [psi,fi,wght,bwght,A,f,fic,eta,beta,fix,psix,hx,ux,B,G,K]=adjssa(N,n,ist,sigma,u,H,bxc,AA,rhoig,dx)
% solve for the adjoint of the SSA equation
% Input
% N=number of inner points, boundaries at 0 and N+1
% n in Glen's flow law
% m in Wertman's friction law
% u velocity from forward eq
% H thickness from forward eq
% betain coeff in front of u in friction law (from forward solve)
% bxc constant inclination of bottom surface
% AA constant in viscosty (from forward solve)
% rhoig rho for ice times g
% h step length
% x0, xN first and last point of interval (not used)
% Output
% psi and fi in adjoint eq
% wght weight in perturbation integral
% A system matrix
% fic coeff multiplying fi_x in 1st eq
    C=2.4126e4;
    A2=2/AA^(1/n);
              m=1/3;
    %                          m=1;
             mm1=m-1;
    gamc=0;
    gamc1=0;
    gamc2=0;
    B=zeros(N,10);
    G=zeros(N,8);
    K=zeros(N,3);
    ux=zeros(N,1);
    eta=ux;
    beta=ux;
    fic=ux;
    psi=ux;
    fi=ux;
    f=ux;
    bx=ux;
    gam=ux;
    gam1=gam;
    gam2=gam;
    NN=2*N;
    A=zeros(NN,NN);
    D=eye(NN);
    psifi=zeros(NN,1);
    Fh=ux;
    Fu=ux;
    %
    p=(1-n)/n;
    N1=N-1;
    hh=2*dx;
    h2=dx^2;
    % u derivative
    for i=2:N1
        ux(i)=(u(i+1)-u(i-1))/hh;
    end
    ux(1)=ux(2);
    ux(N)=ux(N1);
    % artificial damping
    Ngam=ceil(0.99*N);
    for i=Ngam:N
        gam(i)=gamc*(i-Ngam)/(N-Ngam);
        gam1(i)=gamc1*(i-Ngam)/(N-Ngam);
        gam2(i)=gamc2*(i-Ngam)/(N-Ngam);
    end
    % grounding line index
 %   i=1;
 %   while abs(beta(i))>1e-6
 %        i=i+1;
 %   end
 %   iGL=N;
 %   Fu(iGL)=1/h;
%    Fu(iGL)=H(iGL)/h;
%    Fh(iGL)=u(iGL)/h;
    gaussfact=1/sqrt(2*pi*sigma^2);
    for i=1:N
        xi(i)=(i-ist)*dx;
        Fu(i)=gaussfact*exp(-xi(i)^2/(2*sigma^2));
    end
    Fint=dx*sum(Fu)
    Fu=Fu./Fint;
    % viscosity, fi coeff in first eq, bottom derivative
    for i=1:N
        eta(i)=A2*abs(ux(i))^p;
        fic(i)=eta(i)*ux(i)-rhoig*H(i);
    end
    for i=1:ist
        beta(i)=C*abs(u(i))^mm1;
    end
    for i=1:N
        bx(i)=bxc;
    end
    % Glen's coefficient in adjoint viscosity
          nn=n;
          %       nn=1;
    % first point
    i=1;
    j=2*i-1;
    k=2*i;
    jp=2*i+1;
    kp=2*i+2;
    % 1st eq
    A(j,j)=-u(i)/dx;
    A(j,jp)=u(i)/dx;
    A(j,k)=-fic(i)/dx+bx(i)*rhoig;
    A(j,kp)=fic(i)/dx;
    f(j)=Fh(i);
    % 2nd eq
    A(k,j)=H(i)/dx;
    A(k,jp)=-H(i)/dx;
    etaHp=(eta(i+1)*H(i+1)+eta(i)*H(i))/(2*nn);
    A(k,k)=-2*etaHp/h2-beta(i)*m;
    A(k,kp)=etaHp/h2;
        B(i,1)=H(i)/dx;
        B(i,2)=-H(i)/dx;
        B(i,4)=-2*etaHp/h2;
        B(i,5)=-beta(i)*m;
        B(i,3)=0;
        B(i,6)=etaHp/h2;
        B(i,7)=-Fu(i);
    f(k)=-Fu(i);
%    A(k,k)=1;
    f(k)=0;
    % inner points
    for i=2:N1
        j=2*i-1;
        k=2*i;
        jp=2*i+1;
        jm=2*i-3;
        kp=2*i+2;
        km=2*i-2;
        % 1st eq
        A(j,jm)=gam1(i)/dx;
        A(j,j)=-u(i)/dx-2*gam1(i)/dx;
        A(j,jp)=u(i)/dx+gam1(i)/dx;
        A(j,km)=gam2(i)/dx;
        A(j,k)=-fic(i)/dx+bx(i)*rhoig-2*gam2(i)/dx;
        A(j,kp)=fic(i)/dx+gam2(i)/dx;
        G(i,1)=-u(i)/dx;
        G(i,2)=u(i)/dx;
        G(i,3)=-fic(i)/dx;
        G(i,4)=fic(i)/dx;
        G(i,5)=rhoig*bx(i);
        f(j)=Fh(i);
        % 2nd eq
        A(k,jm)=gam(i)/dx;
        A(k,j)=H(i)/dx-2*gam(i)/dx;
        A(k,jp)=-H(i)/dx+gam(i)/dx;
        etaHm=(eta(i-1)*H(i-1)+eta(i)*H(i))/(2*nn);
        etaHp=(eta(i+1)*H(i+1)+eta(i)*H(i))/(2*nn);
        A(k,k)=-(etaHm+etaHp)/h2-beta(i)*m;
        A(k,km)=etaHm/h2;
        A(k,kp)=etaHp/h2;
        B(i,1)=H(i)/dx;
        B(i,2)=-H(i)/dx;
        B(i,4)=-(etaHm+etaHp)/h2;
        B(i,5)=-beta(i)*m;
        B(i,3)=etaHm/h2;
        B(i,6)=etaHp/h2;
        B(i,7)=-Fu(i);
        f(k)=-Fu(i);
    end
    % last point
    i=N;
    j=2*i-1;
    k=2*i;
    jm=2*i-3;
    km=2*i-2;
    % 1st eq
    A(j,j)=-u(i)/dx-2*gam1(i)/dx;
    A(j,km)=gam2(i)/dx;
    A(j,k)=-fic(i)/dx+bx(i)*rhoig-2*gam2(i)/dx;
    f(j)=Fh(i);
    % 2nd eq
    A(k,jm)=gam(i)/dx;
    A(k,j)=H(i)/dx-2*gam(i)/dx;
    etaHm=(eta(i-1)*H(i-1)+eta(i)*H(i))/(2*n);
    A(k,k)=-2*etaHm/h2-beta(i)*m;
    A(k,km)=etaHm/h2;
        B(i,1)=H(i)/dx;
        B(i,2)=0;
        B(i,4)=-2*etaHm/h2;
        B(i,5)=-beta(i)*m;
        B(i,3)=etaHm/h2;
        B(i,6)=0;
        B(i,7)=-Fu(i);
    f(k)=-Fu(i);
    %    HGL=H(iGL);
    %    viGL=eta(iGL)*HGL/nn;
    %    HfiGL=HGL*fic(iGL)/u(iGL);
    %    A(j,j)=-HGL;
    %    A(j,k)=viGL/h;
    %    A(j,km)=-viGL/h;
    %    f(j)=1;
    %    A(k,k)=-viGL/h-HfiGL;
    %    A(k,km)=viGL/h;
    %    f(k)=-1;
%    A(k,k)=1;
%    A(k,k)=fic(iGL);
%    A(k,j)=u(iGL);
%    f(k)=0;
    % rescale A
    Ascal=1e6;
    for i=1:2:NN-1
        D(i,i)=Ascal*D(i,i);
    end
    A=A*D;
    % solve for psi, fi
    psifi=A\f;
    % extract psi, fi
    for i=1:N
        j=2*i-1;
        k=2*i;
        psi(i)=psifi(j);
        fi(i)=psifi(k);
    end
    psi=psi*Ascal;
    for i=2:N1
        fix(i)=(fi(i+1)-fi(i-1))/hh;
        hx(i)=(H(i+1)-H(i-1))/hh+bx(i);
        psix(i)=(psi(i+1)-psi(i))/dx;
    end
    fix(1)=fix(2);
    psix(1)=psix(2);
    hx(1)=hx(2);
    fix(N)=fix(N1);
    psix(N)=psix(N1);
    hx(N)=hx(N1);
    for i=2:N1
        B(i,8)=B(i,3)*fi(i-1)+B(i,4)*fi(i)+B(i,6)*fi(i+1);
        B(i,9)=B(i,5)*fi(i);
        B(i,10)=B(i,1)*psi(i)+B(i,2)*psi(i+1);
        G(i,6)=G(i,2)*psi(i+1)+G(i,1)*psi(i);
        G(i,7)=G(i,4)*fi(i+1)+G(i,3)*fi(i);
        G(i,8)=G(i,5)*fi(i);
        K(i,1)=nn*(B(i,3)*u(i-1)+B(i,4)*u(i)+B(i,6)*u(i+1));
        K(i,2)=-beta(i)*u(i);
        K(i,3)=-rhoig*H(i)*hx(i);
    end
    % weight in perturbed integral
    for i=1:N
        wght(i)=-fi(i)*u(i)^m;
        bwght(i)=psix(i)*u(i)+fix(i)*eta(i)*ux(i)+rhoig*fi(i)*hx(i);
    end
    
    
         
   