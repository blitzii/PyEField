%Testing

lengthscale = 0.1e-6;

%Make the matrix

phi_in = zeros(51,2201);

% Bottom zero
phi_in(end,:) = 0;
% Top electrode 10
phi_in(1,:) = 0;
% Bottom electrode (patterned) -10

% Original (manual/static) pattern

%pattern = zeros(1,101);
%pattern(end-5:end) = -10;

% New (adjustable) pattern




pattern = zeros(1,441);
pattern(101:121) = -10;
pattern(end-120:end-100) = 10;

repeats = round(length(phi_in)/length(pattern));

rar = repmat(pattern,1,repeats);
rar = rar(1:end-(repeats-1));


%phi_in(end-5,:) = rar;
%phi_in(end-4,:) = rar;
%phi_in(end-3,:) = rar;
%phi_in(end-2,:) = rar;
%phi_in(end-1,:) = rar;
phi_in(end,:) = rar;

rar2 = abs(rar./-10);


mask = zeros(51,2201);
mask(end,:) = rar2;
%mask(end-5,:) = rar2;
%mask(end-4,:) = rar2;
%mask(end-3,:) = rar2;
%mask(end-2,:) = rar2;
%mask(end-1,:) = rar2;


mask(1,:) = 1;
mask(end,:) = 1;

permgrid = ones(50,2200);


permgrid = permgrid*2; % Dodecane in bulk
%permgrid = permgrid*80; % Water in bulk
permgrid(1,:) = 6; %Glass at top
permgrid(end,:) = 6; %Glass at bottom

[phi_out] = FieldSolverMG(phi_in,mask,permgrid,1e-6,4);





% Electric field calculation
[Ex,Ez] = phi2E(phi_out{1});

i = 1;
j = 1;

[sizey,sizex] = size(phi_in);

% This correction due to how size is retrieved
sizey = sizey -1;

while i < sizex + 1
   
    E(j,i) = (Ex(j,i).^2 + Ez(j,i).^2).^0.5;
    j = 1;
    
    while j < sizey + 1
        
        E(j,i) = (Ex(j,i).^2 + Ez(j,i).^2).^0.5;
        j = j + 1;
        
    end
    
    i = i + 1;
    
end


% Set E = 0 where mask is
 antimask = mask < 1;
 E = E .* antimask;

%E = E * (1e-6/lengthscale);

% Plotting!
subplot(2,5,1)
figure(1)
contourf(phi_in)
title('Voltage in')

subplot(2,5,2)
contourf(mask)
title('Masked regions')

subplot(2,5,3)
contourf(permgrid)
title('Permitivitty grid')

subplot(2,5,4)
contourf(phi_out{1})
title('Voltage Profile')

subplot(2,5,5)
contourf(E,10)
title('Electric Field Strength')

