function [phiout,deltaphi] = FieldRelaxTol(phi,mask,epsilin,tolerance,quiet,residual)
% FieldRelaxTol: A function that solves for the 2D field around specified
%                electrode patterns until a tolerance limit has been 
%                reached. The potential at the upper and lower surfaces 
%                must be constrained, and side boundary conditions are 
%                reflective.
%
% 
% Usage example:
%  >> [phiout, deltaphi] = FieldRelaxTol(phi,mask,epsilin,tolerance,quiet,residual);
%
% Inputs:
% phi:        A matrix specifying the starting potential
% mask:       The electrode mask: a matrix of the same dimensions as phi in
%             which the elements are either 1 or 0. At the positions where
%             the mask is 1, the potential is fixed at that specified in
%             phi. Where mask is 0, the potential is free to change (make
%             mask 1 where the electrodes are)
% epsilin:    A matrix specifying the relative permittivity of the system.
%             This should be one element shorter than phi in both
%             dimensions.
% tolerance:  The tolerance level. The program will iterate and calculate
%             the potential until the maximum change in the potential in
%             one iteration is less than this value.
% quiet:      If 0, the routine displays the maximum change in phi during 
%             the final iteration (set to zero if note specified).
% residual:   A matrix the same size as phi that specifies the residual
%             from a previous calculation (set to zero if not specified).
%
% Outputs:
% phiout:     The calculated potential
% deltaphi:   The change in phi calculated for the final iteration.
%
%     Adjusted for circular symmetry (I think!).
%
%     Paul Brimicombe, 2010.

%time = cputime;
[zsize, xsize] = size(phi);
phiout = phi;

if nargin<6
    residual = zeros(zsize,xsize);
end

if nargin<5
    quiet = 1;
end

if mean(mask(1,:)) && mean(mask(end,:))
    interdig = 0;
else
    interdig = 1;
end

% Define the sizes of the arrays required
deps_dz = zeros(zsize,xsize);
deps_dx = zeros(zsize,xsize);
deltaphi = ones(zsize,xsize);

epsil = ones(zsize,xsize);
% Interpolate the input permittivities onto the potential grid
epsil(2:end-1,2:end-1) = (epsilin(2:end,1:end-1)+epsilin(1:end-1,1:end-1)+epsilin(2:end,2:end)+epsilin(1:end-1,2:end))./4;
epsil(1,:) = epsil(2,:);
epsil(end,:) = epsil(end-1,:);
epsil(:,1) = epsil(:,2);
epsil(:,end) = epsil(:,end-1);

% Differentiate epsilon wrt z
deps_dz(2:end-1,2:end-1) = (epsilin(2:end,1:end-1)+epsilin(2:end,2:end) - epsilin(1:end-1,1:end-1) - epsilin(1:end-1,2:end))./2;
deps_dz(:,1) = deps_dz(:,2);
deps_dz(:,end)=deps_dz(:,end-1);

% Differentiate epsilon wrt x.
deps_dx(2:end-1,2:end-1) = (epsilin(1:end-1,2:end)+epsilin(2:end,2:end) - epsilin(1:end-1,1:end-1) - epsilin(2:end,1:end-1))./2;
deps_dx(:,1) = 0;
deps_dx(:,end) = 0;
    
% Pre-calculate the epsilon factor (saves us doing it every time!).
epsfact = 1./(2.*epsil);

% The amount of SOR used: decrease if unstable!
sigma = 1.3;
%sigma = 0.2;

iterations=0;
% Loop and update the potential, satisfying the boundary conditions
while max(max(abs(deltaphi(2:end-1,2:end-1)))) > tolerance
    iterations=iterations+1;
     % Use red and black ordering
     % First set of squares
    deltaphi(2:2:end-1,2:2:end-1) = ~mask(2:2:end-1,2:2:end-1).*(1/4).*(phiout(3:2:end,2:2:end-1)+phiout(1:2:end-2,2:2:end-1)+phiout(2:2:end-1,3:2:end)+phiout(2:2:end-1,1:2:end-2)-4.*phiout(2:2:end-1,2:2:end-1)+...
        epsfact(2:2:end-1,2:2:end-1).*(deps_dz(2:2:end-1,2:2:end-1).*(phiout(3:2:end,2:2:end-1)-phiout(1:2:end-2,2:2:end-1))+deps_dx(2:2:end-1,2:2:end-1).*(phiout(2:2:end-1,3:2:end) - phiout(2:2:end-1,1:2:end-2)))-...
        residual(2:2:end-1,2:2:end-1));
    deltaphi(3:2:end-1,3:2:end-1) = ~mask(3:2:end-1,3:2:end-1).*(1/4).*(phiout(4:2:end,3:2:end-1)+phiout(2:2:end-2,3:2:end-1)+phiout(3:2:end-1,4:2:end)+phiout(3:2:end-1,2:2:end-2)-4.*phiout(3:2:end-1,3:2:end-1)+...
        epsfact(3:2:end-1,3:2:end-1).*(deps_dz(3:2:end-1,3:2:end-1).*(phiout(4:2:end,3:2:end-1)-phiout(2:2:end-2,3:2:end-1))+deps_dx(3:2:end-1,3:2:end-1).*(phiout(3:2:end-1,4:2:end) - phiout(3:2:end-1,2:2:end-2)))-...
        residual(3:2:end-1,3:2:end-1));
    
    % Find new phi using Successive Over Relaxation (SOR) method
    phiout(2:2:end-1,2:2:end-1) = (phiout(2:2:end-1,2:2:end-1) + deltaphi(2:2:end-1,2:2:end-1)).*sigma+(1-sigma).*phiout(2:2:end-1,2:2:end-1);
    phiout(3:2:end-1,3:2:end-1) = (phiout(3:2:end-1,3:2:end-1) + deltaphi(3:2:end-1,3:2:end-1)).*sigma+(1-sigma).*phiout(3:2:end-1,3:2:end-1);
    
    % Fix the z boundaries to the input voltages
    phiout(1,:) = phi(1,:);
    phiout(end,:) = phi(end,:);
    
     % Fix x-boundaries
    phiout(:,1) = phiout(:,2);
    phiout(:,end) = phiout(:,end-1);
    
     % Second set of squares
    deltaphi(3:2:end-1,2:2:end-1) = ~mask(3:2:end-1,2:2:end-1).*(1/6).*(phiout(4:2:end,2:2:end-1)+phiout(2:2:end-2,2:2:end-1)+phiout(3:2:end-1,3:2:end)+phiout(3:2:end-1,1:2:end-2)-4.*phiout(3:2:end-1,2:2:end-1)+...
        epsfact(3:2:end-1,2:2:end-1).*(deps_dz(3:2:end-1,2:2:end-1).*(phiout(4:2:end,2:2:end-1)-phiout(2:2:end-2,2:2:end-1))+deps_dx(3:2:end-1,2:2:end-1).*(phiout(3:2:end-1,3:2:end) - phiout(3:2:end-1,1:2:end-2)))-...
        residual(3:2:end-1,2:2:end-1));
    deltaphi(2:2:end-1,3:2:end-1) = ~mask(2:2:end-1,3:2:end-1).*(1/6).*(phiout(3:2:end,3:2:end-1)+phiout(1:2:end-2,3:2:end-1)+phiout(2:2:end-1,4:2:end)+phiout(2:2:end-1,2:2:end-2)-4.*phiout(2:2:end-1,3:2:end-1)+...
        epsfact(2:2:end-1,3:2:end-1).*(deps_dz(2:2:end-1,3:2:end-1).*(phiout(3:2:end,3:2:end-1)-phiout(1:2:end-2,3:2:end-1))+deps_dx(2:2:end-1,3:2:end-1).*(phiout(2:2:end-1,4:2:end) - phiout(2:2:end-1,2:2:end-2)))-...
        residual(2:2:end-1,3:2:end-1));
    
    % Find new phi using Successive Over Relaxation (SOR) method
    phiout(3:2:end-1,2:2:end-1) = (phiout(3:2:end-1,2:2:end-1) + deltaphi(3:2:end-1,2:2:end-1)).*sigma+(1-sigma).*phiout(3:2:end-1,2:2:end-1);
    phiout(2:2:end-1,3:2:end-1) = (phiout(2:2:end-1,3:2:end-1) + deltaphi(2:2:end-1,3:2:end-1)).*sigma+(1-sigma).*phiout(2:2:end-1,3:2:end-1);
    
    % Fix the z boundaries to the input voltages
    phiout(1,:) = phi(1,:);
    phiout(end,:) = phi(end,:);

    % Fix x-boundaries
    phiout(:,1) = phiout(:,2);
    phiout(:,end) = phiout(:,end-1);
        
end
if ~quiet; disp(strcat('Error before residual correction: ',num2str(max(max(abs(deltaphi(2:end-1,2:end-1))))))), disp(strcat('Interations: ',num2str(iterations))), end
%cputime-time