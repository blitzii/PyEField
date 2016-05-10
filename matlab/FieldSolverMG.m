function [phi] = FieldSolverMG(phi_input,mask_input,epsilin_input,tolerance,levmax)
% FieldSolverMG: A function that solves for the 2D field around specified
%                electrode patterns using a multigrid technique until a 
%                tolerance limit has been reached. The potential at the 
%                upper and lower surfaces must be constrained, and side 
%                boundary conditions are reflective.
%
% 
% Usage example:
%  >> [phi] = FieldSolverMG(phi_input,mask_input,epsilin_input,tolerance,levmax);
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
% tolerance:     The tolerance level. The program will iterate and
%                calculate the potential until the maximum change in the 
%                potential in one iteration is less than this value.
% levmax:        The number of coarsening steps to use. For each coarsening
%                step the grid size is halved in each direction. Increasing
%                levmax will can decrease the computation time
%                significantly.
%
% Outputs:
% phi:           The calculated potential field as a cell. Each cell
%                element corresponds to a different grid level, phi{1} is
%                the finest grid and phi{levmax} is the coarsest grid.
%
%     Paul Brimicombe, 2010.

% Store the current time so we can work out how long the calculation took.
time = cputime;
% Make sure that levmax is a whole number
levmax = round(levmax);

% Store the input potential, mask and permittivity
phi{1} = phi_input;
phi_in{1} = phi{1};
mask{1} = mask_input;
epsilin{1} = epsilin_input;

% Store the size of the potential field
[zsize{1},xsize{1}] = size(phi{1});

% The number of iterations to use when correcting the residual
ResidIts = 10;

% Coarsen the grid to maximum size
if levmax ~= 1
    for lev=2:levmax
        % Use the grid coarsening routine to decrease the grid size
        [phi{lev} mask{lev} epsilin{lev}] = coarsen_grid(phi{lev-1},mask{lev-1},epsilin{lev-1});
        % Store the size of the grid at each level
        [zsize{lev} xsize{lev}] = size(phi{lev});
        % Store the input values of the potential for use later when we check
        % the boundary conditions are met
        phi_in{lev} = phi{lev};
    end

    % Uncomment this if you want to check what the coarsened input potential
    % fields look like.
    %return

    % Refine repeatedly and relax each time
    for lev = levmax:-1:2
        disp(strcat('Grid level:',num2str(lev),'...'));
        % Relax at level lev
        [phi{lev} deltaphi] = FieldRelaxTol(phi{lev},mask{lev},epsilin{lev},tolerance,0);
        % Relax on the residual at level lev
        error = FieldRelaxIt(zeros(zsize{lev},xsize{lev}),mask{lev},epsilin{lev},ResidIts,0,-1.*deltaphi);
        % Apply residual correction at level lev
        phi{lev} = error + phi{lev};
                
        % Refine the grid to the next level
        phi{lev-1} = refine_grid(phi{lev},zsize{lev-1},xsize{lev-1});
        % Ensure that the electrodes have the correct potential
        phi{lev-1} = ~mask{lev-1}.*phi{lev-1}+mask{lev-1}.*phi_in{lev-1};
    end
end
disp('Grid level:1...');
% Relax at the finest level
[phi{1} deltaphi] = FieldRelaxTol(phi{1},mask{1},epsilin{1},tolerance,0);
error = FieldRelaxIt(zeros(zsize{1},xsize{1}),mask{1},epsilin{1},ResidIts,0,-1.*deltaphi);
phi{1} = error + phi{1};

% Display the time (in seconds) the the calculation took
disp(strcat('Total simulation time:',num2str(cputime-time),'seconds'));

