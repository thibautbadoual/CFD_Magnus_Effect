function pnew = Poiss_p(N,M,dx,dy,f) 

  

%===============================================
%
%  routine solves pressure Poisson equation 
%  (with Neumann boundary conditions) 
%
%===============================================

NX=N-2;
NY=M-2;

%%% LAPLACIAN WITH NEUMANN BCs (pressure laplace operator)

%%% derivee seconde en x
a11=-1;   % Neumann
b11=-1;   % Neumann
C1 = [1 a11 0]; 
C2 = ones(NX-2,1)*[1 -2 1];
C3 = [0 b11 1];
DXX = spdiags([C1; C2; C3],-1:1,NX,NX)'/dx^2;  %'

clear C1;
clear C2;
clear C3;

%%% derivee seconde en y
a11=-1;   % Neumann
b11=-1;   % Neumann
C1 = [1 a11 0]; 
C2 = ones(NY-2,1)*[1 -2 1];
C3 = [0 b11 1];
DYY = spdiags([C1; C2; C3],-1:1,NY,NY)'/dy^2;  %'

clear C1;
clear C2;
clear C3;

LAPLN =  kron(speye(NY),DXX) + kron(DYY,speye(NX));

% Fixes the value at one point to ensure uniqueness
LAPLN(1,2:end)=0.0;




% Enforce the solvability condition with Neumann BCs

mean_val=sum(sum(f(2:N-1,2:M-1)))/(N-2)/(M-2);
f(2:N-1,2:M-1) = f(2:N-1,2:M-1) - mean_val;

%%% Membre de droite
f2D=f(2:N-1,2:M-1);
f1D = reshape(f2D,NX*NY,1);

%%% Resolution
u1D = LAPLN\f1D;

u2D=reshape(u1D,NX,NY);

pnew=zeros(N,M);
pnew(2:N-1,2:M-1)=u2D;

% ghost cell mapping
pnew(1,:)=pnew(2,:);
pnew(end,:)=pnew(end-1,:);
pnew(:,1)=pnew(:,2);
pnew(:,end)=pnew(:,end-1);
