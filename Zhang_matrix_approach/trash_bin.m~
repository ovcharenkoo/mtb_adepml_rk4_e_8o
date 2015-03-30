close all;  %close all extra windows
clc;  %clear console
clear all; %clear all variables
shft=5;

% total number of grid points in each direction of the grid
NX = 100;
NY = 100;

% size of a grid cell
DELTAX = 10.d0;
DELTAY = DELTAX;

nx_vec=1:NX*DELTAX;	%[m]
ny_vec=1:NY*DELTAY;

% P-velocity, S-velocity and density
cp = 3300.d0;
cs =  1905.d0;
density = 2800.d0;

% total number of time steps
% the time step is twice smaller for this fourth-order simulation,
% therefore let us foruble the number of time steps to keep the same total duration
NSTEP =  700;

% time step in seconds
% 8th-order in space and 4th-order in time finite-difference schemes
% are less stable than second-order in space and second-order in time,
% therefore let us divide the time step by 2
DELTAT = 1.d-3;

% display information on the screen from time to time
% the time step is twice smaller for this fourth-order simulation,
% therefore let us foruble the interval in time steps at which we display information
IT_DISPLAY = 10;
% value of PI
PI = 3.141592653589793238462643d0;

%Use explosive source or gaussian?
EXPLOSIVE_SOURCE=true;
% parameters for the explosivesource
fexp = 20.d0;
texp = 1.20d0 / fexp;
factorexp = 1.d4;

% source
ISOURCE = NX - round(NX/3)+shft;
JSOURCE = round(NY / 3) + 1+shft;

xsource = (ISOURCE - 1) * DELTAX;
ysource = (JSOURCE - 1) * DELTAY;
% angle of source force clockwise with respect to vertical (Y) axis
ANGLE_FORCE = 135.d0;

% main arrays
lambda=zeros((NX+4+shft),(NY+4+shft));
mu=zeros((NX+4+shft),(NY+4+shft));
rho=zeros((NX+4+shft),(NY+4+shft));

vx=zeros((NX+4+shft),(NY+4+shft));
vy=zeros((NX+4+shft),(NY+4+shft));
sigmaxx=zeros((NX+4+shft),(NY+4+shft));
sigmayy=zeros((NX+4+shft),(NY+4+shft));
sigmaxy=zeros((NX+4+shft),(NY+4+shft));

% compute the Lame parameters and density
for j = 1:(NY+4+shft)
    for i = 1:(NX+4+shft)
        rho(i,j) = density;
        mu(i,j) = density*cs*cs;
        lambda(i,j) = density*(cp*cp - 2.d0*cs*cs);
    end
end

%Expansion coefficients(Hixon,1997)
am1=-0.30874;
a0=-0.636;
a1=1.2330;
a2=-0.334;
a3=0.04168;

%RK coefficients
alp=zeros(4);
alp(1)=0.d0;
alp(2)=0.5d0;
alp(3)=0.5d0;
alp(4)=1.d0;

beta=zeros(4);
beta(1)=1.d0/6.d0;
beta(2)=1.d0/3.d0;
beta(3)=1.d0/3.d0;
beta(4)=1.d0/6.d0;

%Forward and backward difference operators
WFdksi=@(i,j,f) (am1*f(i-1,j)+a0*f(i,j)+a1*f(i+1,j)+a2*f(i+2,j)+a3*f(i+3,j))/DELTAX;
WBdksi=@(i,j,f) (-am1*f(i+1,j)-a0*f(i,j)-a1*f(i-1,j)-a2*f(i-2,j)-a3*f(i-3,j))/DELTAX;
WFdeta=@(i,j,f) (am1*f(i,j-1)+a0*f(i,j)+a1*f(i,j+1)+a2*f(i,j+2)+a3*f(i,j+3))/DELTAY;
WBdeta=@(i,j,f) (-am1*f(i,j+1)-a0*f(i,j)-a1*f(i,j-1)-a2*f(i,j-2)-a3*f(i,j-3))/DELTAY;

LFF=@(i,j,f) A*WFdksi(i,j,f)+B*WFdeta(i,j,f);
LBB=@(i,j,f) A*WBdksi(i,j,f)+B*WBdeta(i,j,f);
LFB=@(i,j,f) A*WFdksi(i,j,f)+B*WBdeta(i,j,f);
LBF=@(i,j,f) A*WBdksi(i,j,f)+B*WFdeta(i,j,f);

%Generate grids
sphi='-(1.25*pi*x/max(x)+0.25*pi)'; %sting argument for curved interface
[gx,gy, ksi,eta,J, Ji] = func_curv_jacob_pml(NX,NY,0,0,NX*DELTAX, 0, NY*DELTAY,sphi,DELTAX,DELTAY,false);

%FINISH THIS STUFF WITH DERIVATIVES
dksi_dx=0;
dksi_dy=0;

deta_dx=0;
deta_dy=0;

for it=1:NSTEP
    for i=(1+shft):(NX+shft)
        for j=(1+shft):(NY+shft)
            
            dvx=vx(i,j);
            dvy=vy(i,j);
            dsigmaxx=sigmaxx(i,j);
            dsigmayy=sigmayy(i,j);
            dsigmaxy=sigmaxy(i,j);
            
            U=[dvx, dvy, dsigmaxx, dsigmayy, dsigmaxy]';
            lambdav=lambda(i,j);
            muv=mu(i,j);
            rhov=rho(i,j);
            lambda_plus_two_mu=lambdav+2*muv;

            A=[0 0 dksi_dx/rhov 0 dksi_dy/rhov;
               0 0 0 dksi_dy/rhov dksi_dx/rhov;
               lambda_plus_two_mu*dksi_dx lambdav*dksi_dy 0 0 0;
               lambdav*dksi_dx lambda_plus_two_mu*dksi_dy 0 0 0;
               muv*dksi_dy muv*dksi_dx 0 0 0];

            B=[0 0 deta_dx/rhov 0 deta_dy/rhov;
               0 0 0 deta_dy/rhov deta_dx/rhov;
               lambda_plus_two_mu*deta_dx lambdav*deta_dy 0 0 0;
               lambdav*deta_dx lambda_plus_two_mu*deta_dy 0 0 0
               muv*deta_dy muv*deta_dx 0 0 0];
            
            h1=dt*LFF(i,j,U);
            h2=dt*LBB(i,j,U+alp(2)*h1);
            h3=dt*LFF(i,j,U+alp(3)*h2);
            h4=dt*LBB(i,j,U+alp(4)*h3);
            U=U+beta(1)*h1+beta(2)*h2+beta(3)*h3+beta(4)*h4;
            
            % add the source (force vector located at a given grid point)
            t = (double(it-1)) * DELTAT;
            if EXPLOSIVE_SOURCE
                a = pi*pi*fexp*fexp;
                %source position
                sigmaxx(ISOURCE+shft,JSOURCE+shft)=sigmaxx(ISOURCE+shft,JSOURCE+shft)- factorexp * 2.d0*a*(t-texp)*exp(-a*(t-texp)^2);
                sigmayy(ISOURCE+shft,JSOURCE+shft)=sigmayy(ISOURCE+shft,JSOURCE+shft)- factorexp * 2.d0*a*(t-texp)*exp(-a*(t-texp)^2);
            end
        end
        %------------------------------
        %end of Y loop
    end
    %------------------------------
    %end of X loop
    if(mod(it,IT_DISPLAY) == 0 || it == 5) 
        clf;
        v=vy;
        imagesc(nx_vec,ny_vec,v(1+shft:NX+shft,1+shft:NY+shft)'); 
        xlabel('m');
        ylabel('m');
        colorbar();
        set(gca,'YDir','normal');
        drawnow; hold on;
    end
end
%------------------------------
%end of time loop