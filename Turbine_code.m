%Script for creating an optimal wind turbine and calculating stresses on it

clear
clc
close all
set(0,'defaulttextinterpreter','latex')

filename = 'Cl_cd.xlsx';
Cl_cd = xlsread(filename);

angle_of_attack = Cl_cd(:,1);
Cl_100000 = Cl_cd(:,2);
Cd_100000 = Cl_cd(:,3); 
Cl_200000 = Cl_cd(:,4);
Cd_200000 = Cl_cd(:,5);
Cl_50000 = Cl_cd(:,6);
Cd_50000 = Cl_cd(:,7);


%% Definitions (all units are metric)
TSR = 3.5;
R = 0.177; %total radius of the blade
n = 2; %number of blades
airDense = 1.225;
density = 1040; %density of ABS
v=12; %velocity at which we're expecting max Cp
youngsMod = 2.18e9; %of abs
Hub_radius = 0.0225; %radius of the hub
mu = 1.7965*10^(-5); %dynamic viscosity

momOfInertia = 4.98523914e-5; %for chord length of 1

angVel = (TSR*v)/(R+Hub_radius); %v= air velocity, R= radius

Rpm = angVel*(30/pi);

%Define location of section points on x (all units in meters and kg and m/s)
%0 is at the root
points(:,1) = 0:R/50:R;  %in m

cL = 1.35;

%calculate the relative velocity
velocity = (((points(:,1)+Hub_radius) .* angVel).^2 + ((2/3)*v)^2).^0.5; %calculate the relative wind velocity at all points

%% Airfoil definition
aerofoilDefinition = fopen('airfoil_sg6043.txt', 'r');   %change text to load different set of points
fscanf(aerofoilDefinition,'%c %c %c',3);  %read first line + discard (just header titles)
aeroFoilPoints = fscanf(aerofoilDefinition, '%f %f %f', [3,Inf]); %read the rest of the document and save it 
%into a 3xinf matrix with the rows being x,y and z coordinates respectively
aeroFoilPoints(3,:) = []; %discard the z coordinates of the airfoil

%Area Claculations
foilArea = -simpsonInt(1,length(aeroFoilPoints), aeroFoilPoints', @yFunc); %calculate the area for a unit-chord airfoil


%plot aerofoil
figure(1)
plot(aeroFoilPoints(1,:),aeroFoilPoints(2,:));
axis equal

%% Chord calculations
alpha_1 = atan(R./(TSR*points(1:end,1))); %angle between free-stream wind velocity and rotor plane
chord = zeros(length(points),1);
chord(1:end,1) = ((16*pi*points(1:end,1) .* (sin(alpha_1 ./ 3)).^2 )) ./ (n*cL); %implementation of the Schmitz's rule to find optimal chord

% Plot the original distribution

figure(2)
p1 = plot(points(:,1), chord);
ylabel('Chord length [m]','FontSize',16)
xlabel('Blade radius [m]' , 'FontSize',16);
% titles = sprintf('Graph of lift coefficient against angle of attack for');
title({'Original Schmitz chord distribution'}, 'FontSize',17)
% l = legend([p6, p5, p4, p3, p2,p1], {'$y = 0.0827x + 0.1735$', 'Zero cg shift', '$y = 0.0753x + 0.1717$', 'Forward cg shift of 4in', '$y = 0.0762x + 0.1843$', 'Rearward cg shift of 4in'}, 'Location','best', 'FontSize',15);
% set(l, 'Interpreter', 'latex');
% legend('boxoff')
axis equal
axis([0 points(end,1) 0 0.1])
set(gca,'Color','g')
fig = gcf;
fig.InvertHardcopy = 'off';
% saveas(gcf,'Original_Schmitz.png')

%find chords that are too small and increase them achieving a smooth
%transition to the root
pos = find(chord<=0.042);
pos_max = find(chord==max(chord));
pos2 = pos(pos<pos_max);
pos_middle = round((pos2(end)+pos_max)/2);
k = (chord(pos_middle)-0.042)/points(pos_middle,1)^1.5;
% chord(1:pos_middle) = 0.055:(chord(pos_middle)-0.055)/(pos_middle-1):chord(pos_middle);
chord(1:pos_middle) = k*points(1:pos_middle,1).^1.5 + 0.042;

figure(3)
p1 = plot(points(:,1), chord);
ylabel('Chord length [m]','FontSize',16)
xlabel('Blade radius [m]' , 'FontSize',16);
% titles = sprintf('Graph of lift coefficient against angle of attack for');
title({'Modifier Schmitz chord distribution'}, 'FontSize',17)
% l = legend([p6, p5, p4, p3, p2,p1], {'$y = 0.0827x + 0.1735$', 'Zero cg shift', '$y = 0.0753x + 0.1717$', 'Forward cg shift of 4in', '$y = 0.0762x + 0.1843$', 'Rearward cg shift of 4in'}, 'Location','best', 'FontSize',15);
% set(l, 'Interpreter', 'latex');
% legend('boxoff')
axis equal
axis([0 points(end,1) 0 0.1])
set(gca,'Color','g')
fig = gcf;
fig.InvertHardcopy = 'off';
% saveas(gcf,'Modified_Schmitz.png')


%%
%Figure out the Re of each point on the blade and assign AoA,
%cL and cD respectively!!!

Re = (airDense .* velocity .* chord) ./ (mu);


Re_50000 = find(Re<=50000);
Re_50000_2_100000 = find(Re>50000 & Re<=100000);
Re_100000_2_200000 = find(Re>100000 & Re<=200000);
Re_200000 = find(Re>200000);

AoA(Re_50000) = 8.75;
AoA(Re_50000_2_100000) = (1-(Re(Re_50000_2_100000)-50000)/50000) .* 8.75 + ((Re(Re_50000_2_100000)-50000)/50000) .* 7;
AoA(Re_100000_2_200000) = (1-(Re(Re_100000_2_200000)-100000)/100000) .* 7 + ((Re(Re_100000_2_200000)-100000)/100000) .* 5.5;
AoA(Re_200000) = 5.5;
AoA = AoA';

cD(Re_50000) = spline(angle_of_attack,Cd_50000,AoA(Re_50000));
cD(Re_50000_2_100000) = (1-(Re(Re_50000_2_100000)-50000)/50000) .* spline(angle_of_attack,Cd_50000,AoA(Re_50000_2_100000)) + ((Re(Re_50000_2_100000)-50000)/50000) .* spline(angle_of_attack,Cd_100000,AoA(Re_50000_2_100000));
cD(Re_100000_2_200000) = (1-(Re(Re_100000_2_200000)-100000)/100000) .* spline(angle_of_attack,Cd_100000,AoA(Re_100000_2_200000)) + ((Re(Re_100000_2_200000)-100000)/100000) .* spline(angle_of_attack,Cd_200000,AoA(Re_100000_2_200000));
cD(Re_200000) = spline(angle_of_attack,Cd_200000,AoA(Re_200000));
cD = cD';

cL(Re_50000) = spline(angle_of_attack,Cl_50000,AoA(Re_50000));
cL(Re_50000_2_100000) = (1-(Re(Re_50000_2_100000)-50000)/50000) .* spline(angle_of_attack,Cl_50000,AoA(Re_50000_2_100000)) + ((Re(Re_50000_2_100000)-50000)/50000) .* spline(angle_of_attack,Cl_100000,AoA(Re_50000_2_100000));
cL(Re_100000_2_200000) = (1-(Re(Re_100000_2_200000)-100000)/100000) .* spline(angle_of_attack,Cl_100000,AoA(Re_100000_2_200000)) + ((Re(Re_100000_2_200000)-100000)/100000) .* spline(angle_of_attack,Cl_200000,AoA(Re_100000_2_200000));
cL(Re_200000) = spline(angle_of_attack,Cd_200000,AoA(Re_200000));
cL = cL';

twist_angle = atan((2/3)*((R+Hub_radius)./(TSR*(points(1:end,1)+Hub_radius))))-deg2rad(AoA); %Betz twist

%% Twist and mass calculations
points(:,2) = foilArea * (chord).^2;  %scales the unit area and produces actual area of the airfoil at given chord lengths
Volume = simpsonInt(1, length(points), points, @yFunc); %integrate the function of Chord area vs radius
mass = density*Volume;

%%
%create a new function that incorporates the hub radius into the equation
%allowing for the centrifugal force calculation
r_with_hub = points(:,1) + Hub_radius;

%Centripetal force calculation
centIntFuncHandle = @centIntFunc;
centResult = simpsonInt(1, length(r_with_hub), [r_with_hub, points(:,2)], centIntFuncHandle); %integrate the function of y=radius*area of chord at that radius
centResult = centResult * density * angVel^2;  %total centripetal force at the root

%% centroid calculation for a unit airfoil
[xRoid, yRoid] = centroidCalc(foilArea, aeroFoilPoints);
figure(1)
hold on
scatter(xRoid,yRoid)

%% Calculate the forces per unit span over the blade

lift = cL .* 0.5 .* airDense .* velocity.^2 .* chord;
drag = cD .* 0.5 .* airDense .* velocity.^2 .* chord;

%% Torque calculations
gamma = atan(((1-(1/3))*v)./(angVel*r_with_hub));
F_v = lift .* sin(gamma) - drag .* cos(gamma);
F_w = lift .* cos(gamma) + drag .* sin(gamma);
torque = r_with_hub .* F_v;    %per unit span
torques = [r_with_hub,torque(:,1)];
Total_torque = simpsonInt(1, length(torques), torques, @yFunc);

%% Deflection Calculations
scaledMomOfInertia = momOfInertia * chord .^4;

%vertical Deflec Integrations
yFuncHandle = @yFunc;
forceDistribution = integrateDistribution(points(:,1), F_w, yFuncHandle); %integrate the force/length to get force distribution
Total_force = forceDistribution(end);
forceDistribution = forceDistribution - forceDistribution(end); %correct the force by a constant calculated using cantilever beam theory

r2 = points(:,1)';
r2(length(points(:,1))+1:length(points(:,1))+2) = [points(end,1),0];
w2 = F_w';
w2(length(F_w)+1:length(F_w)+2) = [0,0];
polyin = polyshape({r2},{w2});
[x,y] = centroid(polyin);
% [xRoid, yRoid] = centroidCalc(Total_force, [points(:,1);verticalDistReact]);
% xRoid = -xRoid;
momentDistribution = integrateDistribution(points(:,1), forceDistribution', yFuncHandle);
momentDistribution = momentDistribution + x*Total_force; %see notes
d2vdz2 = -momentDistribution./(scaledMomOfInertia'.*youngsMod);
dvdz = integrateDistribution(points(:,1), d2vdz2', yFuncHandle); %constant is equal to 0
deflection = integrateDistribution(points(:,1), dvdz', yFuncHandle); %constant is equal to 0


figure(4)
p1 = plot(points(:,1), momentDistribution);
ylabel('Bending moment [Nm]','FontSize',16)
xlabel('Blade radius [m]' , 'FontSize',16);
% titles = sprintf('Graph of lift coefficient against angle of attack for');
title({'Bending moment distribution'}, 'FontSize',17)
% axis equal
% axis([0 points(end,1) deflection(end) 0])
set(gca,'Color','g')
fig = gcf;
fig.InvertHardcopy = 'off';
% saveas(gcf,'Bending_moment.png')



figure(5)
p1 = plot(points(:,1), deflection);
ylabel('Blade deflection [m]','FontSize',16)
xlabel('Blade radius [m]' , 'FontSize',16);
% titles = sprintf('Graph of lift coefficient against angle of attack for');
title({'Blade deflection distribution'}, 'FontSize',17)
% l = legend([p6, p5, p4, p3, p2,p1], {'$y = 0.0827x + 0.1735$', 'Zero cg shift', '$y = 0.0753x + 0.1717$', 'Forward cg shift of 4in', '$y = 0.0762x + 0.1843$', 'Rearward cg shift of 4in'}, 'Location','best', 'FontSize',15);
% set(l, 'Interpreter', 'latex');
% legend('boxoff')
% axis equal
axis([0 points(end,1) deflection(end) 0])
set(gca,'Color','g')
fig = gcf;
fig.InvertHardcopy = 'off';
% saveas(gcf,'Deflection.png')




fprintf('Final centrifigual force is %2.2f N \n and stress is %2.2f MPa \n', centResult, max(centResult*1e-6 ./ points(:,2)));

%unit definitions
I_yy = 3.65939038*10^(-3);
Max_x_dist = 0.58422506;
Max_y_dist = 0.0585451461;
I_yy = (momOfInertia*0.7 + 0.3*I_yy);

%scaled definitions
I_yy_scaled = I_yy * chord.^4;
Max_x_dist_scaled = Max_x_dist * chord;
Max_y_dist_scaled = Max_y_dist * chord;

M_w_max = max(momentDistribution);

stress_w = (momentDistribution' .* Max_y_dist_scaled)./ scaledMomOfInertia;



%% vertical Deflec Integrations
yFuncHandle = @yFunc;
force_in_v = integrateDistribution(points(:,1), F_v, yFuncHandle); %integrate the force/length to get force distribution
Total_force_v = force_in_v(end);
force_in_v = force_in_v - force_in_v(end); %correct the force by a constant calculated using cantilever beam theory

r3 = points(:,1)';
r3(length(points(:,1))+1:length(points(:,1))+2) = [points(end,1),0];
w3 = F_v';
w3(length(F_v)+1:length(F_v)+2) = [0,0];
polyin = polyshape({r3},{w3});
[x2,y2] = centroid(polyin);

M_v = integrateDistribution(points(:,1), force_in_v', yFuncHandle);
M_v = M_v + x2*Total_force_v; %see notes

M_v_max = max(M_v);

stress_v = (M_v' .* Max_x_dist_scaled)./ I_yy_scaled;

figure(6)
p1 = plot(points(:,1), stress_v);
ylabel('Normal stress [N/m$^2$]','FontSize',16)
xlabel('Blade radius [m]' , 'FontSize',16);
% titles = sprintf('Graph of lift coefficient against angle of attack for');
title({'Normal stress distribution'}, 'FontSize',17)
% axis equal
% axis([0 points(end,1) deflection(end) 0])
set(gca,'Color','g')
fig = gcf;
fig.InvertHardcopy = 'off';
% saveas(gcf,'Stress.png')

fprintf('Max stress in v and w directions are %2.2f MPa and %2.2f MPa respectively \n', max(stress_v)/(10^6), max(stress_w)/(10^6));
fprintf('\t Total power generated is %2.2f W and the Cp is %2.2f \n', Total_torque*angVel, (Total_torque*angVel)/(0.5*airDense*pi*((r_with_hub(end))^2)*12^3));
function vol = centIntFunc(points)
    vol = points(:,1) .* points(:,2);
end

