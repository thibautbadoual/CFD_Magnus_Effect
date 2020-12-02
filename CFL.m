function vcfl = CFL(u,v)
% computes the CFL condition

global dx;
global dy;

umax = max(max(max(u)),0.01);
vmax = max(max(max(v)),0.01);

vcfl = 0.8*min(dx/umax,dy/vmax);

end