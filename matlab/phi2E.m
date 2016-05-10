function [Ex,Ez] = phi2E(phi)
% Calculate the electric field from a 2D voltage map.
[zsize, xsize] = size(phi);

[Ex,Ez] = gradient(phi,1,1);
Ex = -Ex;
Ez = -Ez;

for i = 2:zsize-1
    Ex(i,1) = (phi(i,xsize)-phi(i,2))./2;
    Ex(i,xsize) = (phi(i,xsize-1)-phi(i,1))./2;
end
