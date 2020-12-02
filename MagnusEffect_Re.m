%==========================================
%
%       Navier-Stokes flow opened domain
%
%==========================================

clear all
close all

% Global variables
global M;
global N;
global dx;
global dy;
global uadv;
global vadv;
global gradphix;
global gradphiy;
global vort;

%Propriété du fluide
L=1;
U=2;
L_Re = [1., 1.e2, 1.e4, 1.e6];
rho = 1.293e-3;

%Propriété de la balle en rotation
R = 0.5;
Omega = 2*pi/5;
m = 10;

niter = 0;
nitermax = 50;

% Domain size
Long = 20*L;
Larg = 10*L;

% Grid size
M = 300;
N = 150;

dx = Long/(M-2);
dy = Larg/(N-2);

% Gridding
x = linspace(dx/2,Long-dx/2,M-2); 
y = linspace(dy/2,Larg-dy/2,N-2); 
[xx,yy] = meshgrid(x,y); 

L_X1 = zeros(nitermax,1);
L_Y1 = zeros(nitermax,1);
L_X2 = zeros(nitermax,1);
L_Y2 = zeros(nitermax,1);
L_X3 = zeros(nitermax,1);
L_Y3 = zeros(nitermax,1);
L_X4 = zeros(nitermax,1);
L_Y4 = zeros(nitermax,1);

for ij=1:4
   
  Re = L_Re(ij)
  niter = 0;
  
  CX = 4*M/5;
  CY = N/2;

  Xi = CX;
  Yi = CY;
  V_S = [-5 0];

  % Initializations
  u = zeros(N,M);
  v = zeros(N,M);

  vort = zeros(N,M);

  uadv = zeros(N,M);
  vadv = zeros(N,M);

  ustar = zeros(N,M);
  vstar = zeros(N,M);

  divstar = zeros(N,M);

  phi = zeros(N,M);
  gradphix = zeros(N,M);
  gradphiy = zeros(N,M);

  L_Cx = zeros(nitermax,1);
  L_Cm = zeros(nitermax,1);
  % ITERATIONS
  while ((niter<nitermax))
    if((mod(niter,100)==0)&&(niter>1))
        fprintf('Iteration: %d', niter);
        fprintf('\n');
    end
    
    % Define the obstacle
      cx = (Xi)*dx;
      cy = (Yi)*dy;
      % Mask for the code (0 in the solid, 1 outside)
      Solid=ones(N,M);
      Solid(2:N-1,2:M-1)=0.5+0.5*sign(sqrt((xx-cx).^2+(yy-cy).^2)-R);
      % Mask for plotting (NaN in the solid, 1 outside)
      Solid2=Solid;
      Solid2(Solid==0)=NaN;
      
      S_theta_sin= zeros(N,M);
      S_theta_cos= ones(N,M)*pi/2;
      comp = 0;
      for k=2:M-1
        for i=2:N-1
          if Solid(i,k) == 1
            if (Solid(i-1,k) == 0)|| (Solid(i+1,k) == 0) || (Solid(i,k-1) == 0) || (Solid(i,k+1) == 0)
              S_theta_sin(i,k) = atan2d(i-Yi, k-Xi)*pi/180;
              S_theta_cos(i,k) = atan2d(i-Yi, k-Xi)*pi/180;
              comp = comp + 1;
            end
          end
        end
      end
    
    
      %surf(xx,yy,S_theta_sin(2:N-1,2:M-1))
      
      dt = CFL(u,v);
      % Advection step
      SemiLagUV(u,v,dt);
      % Diffusion step (explicit)
      ustar = dt/Re*Laplacien(u)+uadv;
      vstar = dt/Re*Laplacien(v)+vadv;
      
      % Boundary conditions
      % LEFT
      ustar(1:N,1)=-ustar(1:N,2) + 2.0*U;
      vstar(1:N,1)=-vstar(1:N,2);
      % TOP
      ustar(N,1:M)=ustar(N-1,1:M);
      vstar(N,1:M)=-vstar(N,1:M);
      % BOTTOM
      ustar(1,1:M)=ustar(2,1:M);
      vstar(1,1:M)=-vstar(2,1:M);
      % RIGHT
      ustar(1:N,M)=ustar(1:N,M-2); 
      vstar(1:N,M)=vstar(1:N,M-2); 

      % Immerged solid
      ustar = ((ustar + R * Omega * sin(S_theta_sin)).* Solid);
      vstar = ((vstar - R * Omega * cos(S_theta_cos)).* Solid);
      
      % Compute divergence
      divstar=div(ustar,vstar);
      %Resolve \Delta phi = div ustar
      phi = Poiss(divstar,N,M);
      
      % Boundary conditions
      % LEFT
      phi(2:N-1,1)=phi(2:N-1,2);
      % RIGHT
      phi(2:N-1,M)=-phi(2:N-1,M-1);
      % TOP
      phi(N,2:M)=phi(N-1,2:M);
      % BOTTOM
      phi(1,2:M)=phi(2,2:M);
      
      % Projection step
      grad(phi); 
      u = ustar + gradphix;
      v = vstar + gradphiy;     
      
      d = div(u,v);
      pnew = Poiss_p(N,M,dx,dy,d/dt);
      p = pnew;
      
      %Force de pression
      Fl = - R*2*pi/comp*[sum(sum(p.*cos(S_theta_cos))) sum(sum(p.*sin(S_theta_sin)))];
      
      %Force de trainee et de Magnus théorique
      
      V_rel = V_S - [U 0];
      L_Cx(niter+1) = Fl(1)/(0.5*rho*pi*R*R*norm(V_rel)^2);
      L_Cm(niter+1) = Fl(2)/(0.5*rho*pi*R*R*norm(V_rel)^2);
      
      %PFD
      V_S = V_S + dt/m * (Fl);
      
      %Avancee du cylindre
      Xi = Xi + dt*V_S(1);
      Yi = Yi + dt*V_S(2);
      
      
      niter = niter+1;
      if (ij=1) 
        
      L_X1(niter) = Xi;
      L_Y1(niter) = Yi;

      end
    
      if (ij=2)
        
      L_X2(niter) = Xi;
      L_Y2(niter) = Yi;
      end
      if (ij=3)
        
      L_X3(niter) = Xi;
      L_Y3(niter) = Yi;
      end
      if (ij=4)
        
      L_X4(niter) = Xi;
      L_Y4(niter) = Yi;
        
      end

    end
    avg_Cx = mean(L_Cx)
    std_Cx = std(L_Cx)
    avg_Cm = mean(L_Cm)
    std_Cm = std(L_Cm)  
 
end

figure(1)
plot(floor(L_X1)*dx,L_Y1*dy, floor(L_X2)*dx,L_Y2*dy, floor(L_X3)*dx,L_Y3*dy, floor(L_X4)*dx,L_Y4*dy);
title('Trajectoire pour différents nombres de Reynolds')
legend('Re = 1','Re = 1e2','Re = 1e4','Re = 1e6');
