function [phi, mask, epsilin] = GridCheck(phi_input,mask_input,epsilin_input,levmax)
% GridCheck: A function that allows you to check whether the grid
%            coarsening process retains the required boundary conditions
%            from the input field
% 
% Usage example:
%  >> [phi, mask, epsilin] = GridCheck(phi_input,mask_input,epsilin_input,levmax);
%
% Inputs:
% phi_input:     A matrix specifying the starting potential
% mask_input:    The electrode mask: a matrix of the same dimensions as phi
%                in which the elements are either 1 or 0. At the positions 
%                where the mask is 1, the potential is fixed at that
%                specified in phi. Where mask is 0, the potential is free
%                to change (make mask 1 where the electrodes are)
% epsilin_input: A matrix specifying the relative permittivity of the
%                system. This should be one element shorter than phi in
%                both dimensions.
% levmax:        The number of coarsening steps to use. For each coarsening
%                step the grid size is halved in each direction.
%
% Outputs:
% phi:           The coarsened field information as a cell: phi{1}
%                corresponds to the finest grid, and phi{levmax} the
%                coarsest grid
% mask:          The coarsened mask information as a cell.
% epsilin:       The coarsened permittivity information as a cell
%     Paul Brimicombe, 2010.

% Make sure that levmax is a whole number
levmax = round(levmax);

% Store the input potential, mask and permittivity
phi{1} = phi_input;
phi_in{1} = phi{1};
mask{1} = mask_input;
epsilin{1} = epsilin_input;

% Store the size of the potential field
[zsize{1},xsize{1}] = size(phi{1});

% Coarsen the grid to maximum size
if levmax ~= 1
    for lev=2:levmax
        % Use the grid coarsening routine to decrease the grid size
        [phi{lev} mask{lev} epsilin{lev}] = coarsen_grid(phi{lev-1},mask{lev-1},epsilin{lev-1});
    end
end
