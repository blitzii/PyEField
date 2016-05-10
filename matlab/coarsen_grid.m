function [phi_2h, mask_2h, epsilin_2h] = coarsen_grid(phi,mask,epsilin)
% Coarsens the current grid by a factor of two, symmetrically about the
% centre-line in the x-direction, and always starting from the base in the
% z-direction.  This routine is called from FieldSolverMG and GridCheck.

% Store the size of the finer grid
[zsize_h, xsize_h] = size(phi);
% Find the centre-line in the x-direction
cx_h = (xsize_h-1)/2+1;

% Find the furthest points in z and x that we want
zfar_2h = 2*floor((zsize_h-1)/2)+1;
xfar_2h = 2*floor((xsize_h-1)/4);

% Pre-allocate the array sizes for speed
phi_2h = zeros(size(1:2:zfar_2h,2),size(cx_h+((-xfar_2h):2:xfar_2h),2));
mask_2h = zeros(size(1:2:zfar_2h,2),size(cx_h+((-xfar_2h):2:xfar_2h),2));

% Copy across the phi and mask data from the finer grid to the coarser grid
% Start at the first line:
phi_2h(1,:) = phi(1,cx_h+((-xfar_2h):2:xfar_2h));
mask_2h(1,:) = mask(1,cx_h+((-xfar_2h):2:xfar_2h));
j=0;
% Now do the rest of the lines in turn
for i=2:2:zfar_2h
    j=j+1;
    phi_2h(j,:) = phi(i,cx_h+((-xfar_2h):2:xfar_2h));
    mask_2h(j,:) = mask(i,cx_h+((-xfar_2h):2:xfar_2h));
    % Check whether we have missed any important information (i.e. an
    % electrode), and any required information to the new grid
    if find(mask(i-1,cx_h+((-xfar_2h):2:xfar_2h)))
        mask_2h(j,:) = mask(i-1,cx_h+((-xfar_2h):2:xfar_2h));
        phi_2h(j,:) = phi(i-1,cx_h+((-xfar_2h):2:xfar_2h));
    end
end
% Ensure that the boundary conditions are retained at the top and bottom.
phi_2h(end,:) = phi(end,cx_h+((-xfar_2h):2:xfar_2h));
mask_2h(end,:) = mask(end,cx_h+((-xfar_2h):2:xfar_2h));

% Pre-allocate the permittivity array
epsil = ones(zsize_h,xsize_h);
% Interpolate the input permittivities onto the potential grid at the finer
% grid level
epsil(2:end-1,2:end-1) = (epsilin(2:end,1:end-1)+epsilin(1:end-1,1:end-1)+epsilin(2:end,2:end)+epsilin(1:end-1,2:end))./4;
epsil(1,:) = epsil(2,:);
epsil(end,:) = epsil(end-1,:);
epsil(:,1) = epsil(:,2);
epsil(:,end) = epsil(:,end-1);

% Copy the required values onto the coarser grid
epsilin_2h = epsil(2:2:(zfar_2h-1),cx_h+((-xfar_2h+1):2:(xfar_2h-1)));
