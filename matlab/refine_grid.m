function [phi_h] = refine_grid(phi_2h,zsize_h,xsize_h)
% Refines the potential onto a grid that is twice as fine by interpolation.
% This routine is called from FieldSolverMG.

% Interpolate the potential to a finer grid
phi_h = interp2(phi_2h);

% Store the size of the new grid
[zsize_new, xsize_new] = size(phi_h);

% Check whether it's the same size as it should be...
if zsize_h ~= zsize_new
    % If not, add another element to fix it.
    phi_h(end+1,:) = phi_h(end,:);
end
if xsize_h ~= xsize_new
    % If not, add some elements to fix it.
    phi_h(:,2:end+1) = phi_h;
    phi_h(:,1) = phi_h(:,2);
    phi_h(:,end+1) = phi_h(:,end);
end