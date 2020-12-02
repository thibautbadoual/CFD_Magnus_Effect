function rst = Laplacien(chp)
% Calcule le laplacien scalaire rst(i,j) du champ scalaire chp(i,j)

global N;
global M;
global dx;
global dy;

rst = zeros(N,M);
coefx = 1/dx/dx;
coefy = 1/dy/dy;
coef0 = 2*(coefx + coefy); 

rst(2:N-1,2:M-1) = (chp(2:N-1,3:M) + chp(2:N-1,1:M-2))*coefx + ... 
                      (chp(3:N,2:M-1) + chp(1:N-2,2:M-1))*coefy-chp(2:N-1,2:M-1)*coef0; 
                  
end