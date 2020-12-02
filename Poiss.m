function [ u2D ] = Poiss( f22D,N,M )
% Solves for Laplacien(phi)=f returns phi

global dx
global dy

% 2nd derivative in x
a11=1;   % 1 or 3 for BC: 1 = Neumann, 3 = Dirichlet
b11=3;   % 1 or 3 for BC: 1 = Neumann, 3 = Dirichlet
C1 = [-1 a11 0]; 
C2 = ones(M-2,1)*[-1 2 -1];
C3 = [0 b11 -1];
DXX = spdiags([C1; C2; C3],-1:1,M,M)'/dx^2;  %'

clear C1;
clear C2;
clear C3;

% 2nd derivative in y
a11=1;   % 1 or 3 for BC: 1 = Neumann, 3 = Dirichlet
b11=1;   % 1 or 3 for BC: 1 = Neumann, 3 = Dirichlet
C1 = [-1 a11 0]; 
C2 = ones(N-2,1)*[-1 2 -1];
C3 = [0 b11 -1];
DYY = spdiags([C1; C2; C3],-1:1,N,N)'/dy^2;

LAPL =  kron(speye(N),DXX) + kron(DYY,speye(M));

% right member term
f1D = reshape(f22D',N*M,1);

% solving
u1D = LAPL\f1D;

u2D=reshape(u1D,M,N)';

end

