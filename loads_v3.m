clc
clear
close all

% file_name = "dircol_sqp_final1_scaled.mat";
% if exist(file_name,"file")
%     load(file_name,"log_x","log_param");
% end

if exist("rho_v_z.mat", "file")
    load("rho_v_z.mat");
end

%   
%   Rocket Properties imported from bottom up sizing
%
%   loadcase.component = [station1 x, axial force1, shear force1y, moment1x, shear force1x, moment1y, station2 x, axial force
%   2, shear force2x, moment2x, shear force2y, moment2y,] - with 1 and 2 being the forward and aft
%   ends, respectivly
%
%
%   loadcase 1 - max q
%   loadcase 2 - max q-alpha
%   loadcase 3 - carriage loads
%
%   max possible wind loading throught the jet stream (~4-8miles, at 275mph)
%   we are launching at ~7.5 miles, which is right in this jet stream 
%
%
%
%
%
% V2 UPDATES:
%    We have 11 nodes, 10 elements - 5 DOF Each
%    DOFs: 1-axial, 2-shear,z 3-moment,y 4-shear,y 5-moment,z

%conversions
psi2Pa = 6894.76;


%initial variables for general sizing
%solid_prop_mass = 1508.0927; ORIGINAL
%liquid_prop_mass = 246.6531; ORIGINAL

solid_prop_mass = 1278.8174; 
liquid_prop_mass = 308.9071; 

Din = 0.58;
vehicle_sizing = bottom_up_sizing(solid_prop_mass,liquid_prop_mass,Din);

%pulling out properties for each component using sizing
%section.part.dims = [length (m), mass (kg), Izz (kgm^2), apx thickness, CG (relative to forward edge of component)]
front_profile = vehicle_sizing.CS(1); %area, m^2
side_profile = vehicle_sizing.CS(2);  %area, m^2
payload_dims = vehicle_sizing.second_stage.payload.dims;
PAF_dims = vehicle_sizing.second_stage.PAF.dims;
PLF_dims = vehicle_sizing.first_stage.PLF.dims;
interstage_dims = vehicle_sizing.first_stage.interstage.dims;
liqcasing_dims = vehicle_sizing.second_stage.liqcasing.dims;
ox_dims = vehicle_sizing.second_stage.ox.dims;
fuel_dims = vehicle_sizing.second_stage.fuel.dims;
ss_engine_dims = vehicle_sizing.second_stage.engine.dims;
ss_wiring_dims = vehicle_sizing.second_stage.wiringgnc.dims;
solid_casing_dims = vehicle_sizing.first_stage.solidcasing.dims;
fs_engine_dims = vehicle_sizing.first_stage.engine.dims;
aftskirt_dims = vehicle_sizing.first_stage.aftskirt.dims;
fs_prop_dims = vehicle_sizing.first_stage.solidprop.dims;
avgnc_dims = vehicle_sizing.first_stage.avgnc.dims;

kg2lb = 2.20462;
m2ft = 3.28084;
    
%max_winds
v_maxwind = 275*0.44704; %mph to m/s

%structural masses and CGs for body loads
PLF_massloc = [PLF_dims(2) + payload_dims(2)+PAF_dims(2),(PLF_dims(2)*PLF_dims(6) +payload_dims(2)*payload_dims(6)+ PAF_dims(2)*PAF_dims(6))/(PLF_dims(2) + payload_dims(2)+PAF_dims(2))]; 
liqcasing_massloc = [liqcasing_dims(2) + ox_dims(2) + fuel_dims(2) + ss_wiring_dims(2), (liqcasing_dims(2)*liqcasing_dims(6) + ox_dims(2)*ox_dims(6) + fuel_dims(2)*fuel_dims(6) + ss_wiring_dims(2)*ss_wiring_dims(6))/(liqcasing_dims(2) + ox_dims(2) + fuel_dims(2) + ss_wiring_dims(2))];
interstage_massloc = [interstage_dims(2) + ss_engine_dims(2), (interstage_dims(2)*interstage_dims(6) + ss_engine_dims(2)*ss_engine_dims(6))/(interstage_dims(2) + ss_engine_dims(2))];
solidcasing_massloc = [solid_casing_dims(2)+fs_prop_dims(2),(solid_casing_dims(2)*solid_casing_dims(6)+fs_prop_dims(2)*fs_prop_dims(6))/(solid_casing_dims(2)+fs_prop_dims(2))];
aftskirt_massloc = [aftskirt_dims(2)+fs_engine_dims(2)+avgnc_dims(2), (aftskirt_dims(2)*aftskirt_dims(6)+fs_engine_dims(2)*fs_engine_dims(6)+avgnc_dims(2)*avgnc_dims(6))/(aftskirt_dims(2)+fs_engine_dims(2)+avgnc_dims(2))];

g = 9.81;

%E, I, A, and L for each part
PLF.E = vehicle_sizing.first_stage.PLF.mat(2);
PLF.I = PLF_dims(3);
PLF.A = pi*((Din/2 + PLF_dims(4))^2-(Din/2)^2);
PLF.L = PLF_dims(1);

liqcasing.E = vehicle_sizing.second_stage.liqcasing.mat(2);
liqcasing.I = liqcasing_dims(3);
liqcasing.A = pi*((Din/2 + liqcasing_dims(4))^2-(Din/2)^2);
liqcasing.L = liqcasing_dims(1);

interstage.E = vehicle_sizing.first_stage.interstage.mat(2);
interstage.I = interstage_dims(3);
interstage.A = pi*((Din/2 + interstage_dims(4))^2-(Din/2)^2);
interstage.L = interstage_dims(1);

solid_casing.E = vehicle_sizing.first_stage.solidcasing.mat(2);
solid_casing.I = solid_casing_dims(3);
solid_casing.A = pi*((Din/2 + solid_casing_dims(4))^2-(Din/2)^2);
solid_casing.L = solid_casing_dims(1);

aftskirt.E = vehicle_sizing.first_stage.aftskirt.mat(2);
aftskirt.I = aftskirt_dims(3);
aftskirt.A = pi*((Din/2 + aftskirt_dims(4))^2-(Din/2)^2);
aftskirt.L = aftskirt_dims(1);

%element lengths
Lels = [PLF.L, (liqcasing.L/2 - vehicle_sizing.dH.liq/2), (liqcasing.L/2- vehicle_sizing.dH.liq/2), interstage.L/2, interstage.L/2, vehicle_sizing.CG - (solid_casing_dims(6) - solid_casing.L/2) - 15/(12*m2ft),15/(12*m2ft) - vehicle_sizing.dH.sol,15/(12*m2ft),solid_casing.L - (30/(12*m2ft)+vehicle_sizing.CG - (solid_casing_dims(6) - solid_casing.L/2) - 15/(12*m2ft))-vehicle_sizing.dH.sol,aftskirt.L + vehicle_sizing.dH.sol];

%element Es
Eels = [PLF.E, liqcasing.E, liqcasing.E, interstage.E, interstage.E, solid_casing.E,solid_casing.E,solid_casing.E,solid_casing.E,aftskirt.E];

%element As
Aels = [PLF.A, liqcasing.A, liqcasing.A, interstage.A, interstage.A, solid_casing.A,solid_casing.A,solid_casing.A,solid_casing.A,aftskirt.A];

%element Is
Iels = [PLF.I, liqcasing.I, liqcasing.I, interstage.I, interstage.I, solid_casing.I,solid_casing.I,solid_casing.I,solid_casing.I,aftskirt.I];

%nodelocs = 
nodeloc = [0,Lels(1),sum(Lels(1:2)),sum(Lels(1:3)),sum(Lels(1:4)),sum(Lels(1:5)),sum(Lels(1:6)),sum(Lels(1:7)),sum(Lels(1:8)),sum(Lels(1:9)),sum(Lels(1:10))];

%%load cases
%PICK YOUR LOAD CASE!! :
load_case = 2;

if load_case == 1 %max q
  %body forces
    nx = 1.6;
    nz = 0.3;

    %distributed load parameters
    alpha = 31;
    beta = 90;
    fpa = 0;

    qpsi = 3.3;
    q = qpsi*6894.76; %Pa
    alt = 13050; %m
    atm = atmo(alt);
    maxqv = 277.226824974998; %m/s
    machnum = maxqv/sqrt(1.4*287*atm.T);
    CD = CD(machnum, alpha);
    shear_aero = CD*q*sind(alpha); %Pa/m^2
    dist_aeroshear = shear_aero*Din; %converts to Pa/m - force per unit length along the beam

    qmax_wind = 0.5*atm.rho*v_maxwind^2;
    axial_wind = -CD*qmax_wind*cosd(beta); %Pa/m^2
    shear_wind = CD*qmax_wind*sind(beta); %Pa/m^2
    dist_windshear = shear_wind*Din; %converts to Pa/m - force per unit length along the beam

    %assign from point loads expected (from
    %engine, nose loads, etc) - put in the syntax as follows: 
    % [node #, force (N or Nm)]
    forces = [1,    -q*front_profile*cosd(alpha)*CD;
              2,    0;
              3,    0;
              4,    0;
              5,    0;
              6,    -(nx + 1)*g*PLF_massloc(1);
              8,    (nz)*g*PLF_massloc(1);
              11,   -(nx + 1)*g*liqcasing_massloc(1)/2;
              13,   (nz)*g*liqcasing_massloc(1)/2;
              16,    -(nx + 1)*g*liqcasing_massloc(1)/2;
              18,   (nz)*g*liqcasing_massloc(1)/2;
              21,   -(nx + 1)*g*interstage_massloc(1)/2;
              23,   (nz)*g*interstage_massloc(1)/2;
              26,   -(nx + 1)*g*interstage_massloc(1)/2;
              28,   (nz)*g*interstage_massloc(1)/2; 
              31,   -(nx + 1)*g*solidcasing_massloc(1)/3;
              33,   (nz)*g*solidcasing_massloc(1)/3;
              %36,   -(nx + 1)*g*solidcasing_massloc(1)/4;
              %38,   (nz)*g*solidcasing_massloc(1)/4;
              41,   -(nx + 1)*g*solidcasing_massloc(1)/3;
              43,   (nz)*g*solidcasing_massloc(1)/3;
              46,   -(nx + 1)*g*solidcasing_massloc(1)/3;
              48,   (nz)*g*solidcasing_massloc(1)/3;
              51,   (nx*g*vehicle_sizing.mass(3) -(nx + 1)*g*aftskirt_massloc(1));
              52,   0;
              53,   0;
              54,   0;
              55,   0];
    %assign specified displacements expected- put in the syntax as follows: 
    % [node #, displacement (m or rad)]
    displacements = [36,0;
                     37,0;
                     38,0;
                     39,0;
                     40,0;];
end

if load_case == 2 %maxqa
  %body forces
    nx = 1.6;
    nz = 0.3;

    %distributed load parameters
    alpha = 34.7;
    beta = 90;
    fpa = 0;

    qapsi = 1.95;
    qalpha = qapsi*6894.76;
    q = qalpha/(alpha*pi/180); %Pa
    alt = 16500; %m
    atm = atmo(alt);
    machnum = 1.6;
    maxqv = machnum*sqrt(1.4*287*atm.T);
    CD = CD(machnum, alpha);
    shear_aero = CD*q*sind(alpha); %Pa/m^2
    dist_aeroshear = shear_aero*Din; %converts to Pa/m - force per unit length along the beam

    qmax_wind = 0.5*atm.rho*v_maxwind^2;
    axial_wind = -CD*qmax_wind*cosd(beta); %Pa/m^2
    shear_wind = CD*qmax_wind*sind(beta); %Pa/m^2
    dist_windshear = shear_wind*Din; %converts to Pa/m - force per unit length along the beam

    %assign from point loads expected (from
    %engine, nose loads, etc) - put in the syntax as follows: 
    % [node #, force (N or Nm)]
    forces = [1,    -q*front_profile*cosd(alpha)*CD;
              2,    0;
              3,    0;
              4,    0;
              5,    0;
              6,    -(nx + 1)*g*PLF_massloc(1);
              8,    (nz)*g*PLF_massloc(1);
              11,   -(nx + 1)*g*liqcasing_massloc(1)/2;
              13,   (nz)*g*liqcasing_massloc(1)/2;
              16,    -(nx + 1)*g*liqcasing_massloc(1)/2;
              18,   (nz)*g*liqcasing_massloc(1)/2;
              21,   -(nx + 1)*g*interstage_massloc(1)/2;
              23,   (nz)*g*interstage_massloc(1)/2;
              26,   -(nx + 1)*g*interstage_massloc(1)/2;
              28,   (nz)*g*interstage_massloc(1)/2; 
              31,   -(nx + 1)*g*solidcasing_massloc(1)/3;
              33,   (nz)*g*solidcasing_massloc(1)/3;
              %36,   -(nx + 1)*g*solidcasing_massloc(1)/4;
              %38,   (nz)*g*solidcasing_massloc(1)/4;
              41,   -(nx + 1)*g*solidcasing_massloc(1)/3;
              43,   (nz)*g*solidcasing_massloc(1)/3;
              46,   -(nx + 1)*g*solidcasing_massloc(1)/3;
              48,   (nz)*g*solidcasing_massloc(1)/3;
              51,   (nx*g*vehicle_sizing.mass(3) -(nx + 1)*g*aftskirt_massloc(1));
              52,   0;
              53,   0;
              54,   0;
              55,   0];
    %assign specified displacements expected- put in the syntax as follows: 
    % [node #, displacement (m or rad)]
    displacements = [36,0;
                     37,0;
                     38,0;
                     39,0;
                     40,0;];
end




if load_case == 3 %carriage just before dropping

    %body forces
    nx = -0.7;
    nz = 1.5;

    %distributed load parameters
    alpha = 0;
    beta = 90;
    fpa = 0;

    alt = 12192; %m
    atm = atmo(alt);
    machnum = 0.85;
    maxqv = machnum*sqrt(1.4*287*atm.T); %m/s
    q = 1/2*atm.rho*maxqv^2;
    CD = CD(machnum, alpha);
    shear_aero = CD*q*sind(alpha); %Pa/m^2
    dist_aeroshear = shear_aero*Din; %converts to Pa/m - force per unit length along the beam


    qmax_wind = 0.5*atm.rho*v_maxwind^2;
    axial_wind = -CD*qmax_wind*cosd(beta); %Pa/m^2
    shear_wind = CD*qmax_wind*sind(beta); %Pa/m^2
    dist_windshear = shear_wind*Din; %converts to Pa/m - force per unit length along the beam

    %assign from point loads expected (from
    %engine, nose loads, etc) - put in the syntax as follows: 
    % [node #, force (N or Nm)]
    forces = [1,    -q*front_profile*cosd(alpha)*CD;
              2,    0;
              3,    0;
              4,    0;
              5,    0;
              6,    -(nx + 1)*g*PLF_massloc(1);
              8,    (nz)*g*PLF_massloc(1);
              11,   -(nx + 1)*g*liqcasing_massloc(1)/2;
              13,   (nz)*g*liqcasing_massloc(1)/2;
              16,    -(nx + 1)*g*liqcasing_massloc(1)/2;
              18,   (nz)*g*liqcasing_massloc(1)/2;
              21,   -(nx + 1)*g*interstage_massloc(1)/2;
              23,   (nz)*g*interstage_massloc(1)/2;
              26,   -(nx + 1)*g*interstage_massloc(1)/2;
              28,   (nz)*g*interstage_massloc(1)/2;
              36,   -(nx + 1)*g*solidcasing_massloc(1)/2;
              38,   (nz)*g*solidcasing_massloc(1)/2;
              46,   -(nx + 1)*g*(solidcasing_massloc(1)+aftskirt_massloc(1));
              48,   (nz)*g*(solidcasing_massloc(1)/2 + aftskirt_massloc(1));
              51,   0;
              52,   0;
              53,   0;
              54,   0;
              55,   0];
    %assign specified displacements expected- put in the syntax as follows: 
    % [node #, displacement (m or rad)]
    displacements = [31,0;
                    32,0;
                    34,0;
                    41,0;
                    42,0;
                    44,0;];
end

if load_case == 4 %carriage in transit

    %body forces
    nx = -0.7;
    nz = 1.5;

    %distributed load parameters
    alpha = 0;
    beta = 90;
    fpa = 0;

    %q = qalpha/(alpha*pi/180); %Pa
    alt = 10668; %m
    atm = atmo(alt);
    maxqv = 750*0.44704;
    machnum = maxqv/sqrt(1.4*287*atm.T); %m/s
    q = 1/2*atm.rho*maxqv^2;
    CD = CD(machnum, alpha);
    shear_aero = CD*q*sind(alpha); %Pa/m^2
    dist_aeroshear = shear_aero*Din; %converts to Pa/m - force per unit length along the beam

    qmax_wind = 0.5*atm.rho*v_maxwind^2;
    axial_wind = -CD*qmax_wind*cosd(beta); %Pa/m^2
    shear_wind = CD*qmax_wind*sind(beta); %Pa/m^2
    dist_windshear = shear_wind*Din; %converts to Pa/m - force per unit length along the beam

    %assign from point loads expected (from
    %engine, nose loads, etc) - put in the syntax as follows: 
    % [node #, force (N or Nm)]
    forces = [1,    -q*front_profile*cosd(alpha)*CD;
              2,    0;
              3,    0;
              4,    0;
              5,    0;
              6,    -(nx + 1)*g*PLF_massloc(1);
              8,    (nz)*g*PLF_massloc(1);
              11,   -(nx + 1)*g*liqcasing_massloc(1)/2;
              13,   (nz)*g*liqcasing_massloc(1)/2;
              16,    -(nx + 1)*g*liqcasing_massloc(1)/2;
              18,   (nz)*g*liqcasing_massloc(1)/2;
              21,   -(nx + 1)*g*interstage_massloc(1)/2;
              23,   (nz)*g*interstage_massloc(1)/2;
              26,   -(nx + 1)*g*interstage_massloc(1)/2;
              28,   (nz)*g*interstage_massloc(1)/2;
              36,   -(nx + 1)*g*solidcasing_massloc(1)/2;
              38,   (nz)*g*solidcasing_massloc(1)/2;
              46,   -(nx + 1)*g*(solidcasing_massloc(1)+aftskirt_massloc(1));
              48,   (nz)*g*(solidcasing_massloc(1)/2 + aftskirt_massloc(1));
              51,   0;
              52,   0;
              53,   0;
              54,   0;
              55,   0];
    %assign specified displacements expected- put in the syntax as follows: 
    % [node #, displacement (m or rad)]
    displacements = [31,0;
                    32,0;
                    34,0;
                    41,0;
                    42,0;
                    44,0;];
end

syms x P L E I A

q = zeros(55,1);
Q = zeros(55,1);

%for each load case, have entered values for the distributed loads, the
%fixed points, etc
conn = zeros(10,10);
for i=1:10
    conn(i,:) = [i*5-4; i*5-3; i*5-2; i*5-1;i*5;i*5+1;i*5+2; i*5+3; i*5+4; i*5+5];
end

Kt = (E*I/L)*[12/L^2, 6/L, -12/L^2, 6/L;
              6/L,    4,   - 6/L,   2;
              -12/L^2,-6/L,12/L^2,  -6/L;
              6/L,    2,   -6/L,    4];
Nt_x = [1-3*x^2/L^2+2*x^3/L^3, x-2*x^2/L + x^3/L^2, 3*x^2/L^2 - 2*x^3/L^3, -x^2/L + x^3/L^2];
Nt_xp = diff(Nt_x,x);
Nt_xpp = diff(Nt_xp,x);
Nt_xppp = diff(Nt_xpp,x);

Ka = (E*A/L)*[1 -1;-1,1];
Na = [1-x/L,x/L];


Q_ar = zeros(length(Q),1);
q_ar = zeros(length(q),1);

%assemble stiffness matrix
K = zeros(55);

% Loop over each element to assemble element stiffness matrices
for i = 1:10
    conn_i = conn(i,:);
    Ka_i = subs(Ka, [E,A,L], [Eels(i),Aels(i),Lels(i)]);
    Kt_i = subs(Kt, [E,I,L], [Eels(i),Iels(i),Lels(i)]);
    %Kt_zi = subs(Kt, [E,I,L], [Eels(i),Iels(i),Lels(i)]);

    K(conn_i(1), conn_i(1)) = K(conn_i(1), conn_i(1)) + Ka_i(1,1);
    K(conn_i(6), conn_i(1)) = K(conn_i(6), conn_i(1)) + Ka_i(2,1);
    K(conn_i(1), conn_i(6)) = K(conn_i(1), conn_i(6)) + Ka_i(1,2);
    K(conn_i(6), conn_i(6)) = K(conn_i(6), conn_i(6)) + Ka_i(2,2);

    for j=1:2
        K(conn_i(2*j),conn_i(2*j)) = K(conn_i(2*j),conn_i(2*j)) + Kt_i(1,1);
        K(conn_i(2*j+1),conn_i(2*j)) = K(conn_i(2*j+1),conn_i(2*j)) + Kt_i(2,1);
        K(conn_i(2*j+5),conn_i(2*j)) = K(conn_i(2*j+5),conn_i(2*j)) + Kt_i(3,1);
        K(conn_i(2*j+6),conn_i(2*j)) = K(conn_i(2*j+6),conn_i(2*j)) + Kt_i(4,1);

        K(conn_i(2*j),conn_i(2*j+1)) = K(conn_i(2*j),conn_i(2*j+1)) + Kt_i(1,2);
        K(conn_i(2*j+1),conn_i(2*j+1)) = K(conn_i(2*j+1),conn_i(2*j+1)) + Kt_i(2,2);
        K(conn_i(2*j+5),conn_i(2*j+1)) = K(conn_i(2*j+5),conn_i(2*j+1)) + Kt_i(3,2);
        K(conn_i(2*j+6),conn_i(2*j+1)) = K(conn_i(2*j+6),conn_i(2*j+1)) + Kt_i(4,2);

        K(conn_i(2*j),conn_i(2*j+5)) = K(conn_i(2*j),conn_i(2*j+5)) + Kt_i(1,3);
        K(conn_i(2*j+1),conn_i(2*j+5)) = K(conn_i(2*j+1),conn_i(2*j+5)) + Kt_i(2,3);
        K(conn_i(2*j+5),conn_i(2*j+5)) = K(conn_i(2*j+5),conn_i(2*j+5)) + Kt_i(3,3);
        K(conn_i(2*j+6),conn_i(2*j+5)) = K(conn_i(2*j+6),conn_i(2*j+5)) + Kt_i(4,3);

        K(conn_i(2*j),conn_i(2*j+6)) = K(conn_i(2*j),conn_i(2*j+6)) + Kt_i(1,4);
        K(conn_i(2*j+1),conn_i(2*j+6)) = K(conn_i(2*j+1),conn_i(2*j+6)) + Kt_i(2,4);
        K(conn_i(2*j+5),conn_i(2*j+6)) = K(conn_i(2*j+5),conn_i(2*j+6)) + Kt_i(3,4);
        K(conn_i(2*j+6),conn_i(2*j+6)) = K(conn_i(2*j+6),conn_i(2*j+6)) + Kt_i(4,4);
    end
end

% Display the global stiffness matrix
%disp(K);

%assemble loads array
for i = 1:10
    load_dist_aero = shearmoment(dist_aeroshear,Nt_x,Lels(i),nodeloc(i));
    load_dist_wind = shearmoment(dist_windshear,Nt_x,Lels(i),nodeloc(i));

    Q_ar(5*i-3)=Q_ar(5*i-3) + load_dist_aero(1);
    Q_ar(5*i-2)=Q_ar(5*i-2) + load_dist_aero(2);
    Q_ar(5*i+2)=Q_ar(5*i+2) + load_dist_aero(3);
    Q_ar(5*i+3)=Q_ar(5*i+3) + load_dist_aero(4);
    Q_ar(5*i-1)=Q_ar(5*i-1) + load_dist_wind(1);
    Q_ar(5*i)  =Q_ar(5*i) + load_dist_wind(2);
    Q_ar(5*i+4)=Q_ar(5*i+4) + load_dist_wind(3);
    Q_ar(5*i+5)=Q_ar(5*i+5) + load_dist_wind(4);
end

list = 1:length(Q);

arr_unknown_Q = [];
arr_unknown_q = [];

for i = 1:length(q_ar)
    if ismember(i,displacements(:,1))
        q_ar(i) = q_ar(i) + displacements(find(displacements(:,1)==i),2);
    else
        q_ar(i) = NaN;
        arr_unknown_q = [arr_unknown_q,i];
    end
end


for i = 1:length(Q_ar)
    if ismember(i, forces(:,1))
        if forces(find(forces(:,1)==i),2) == 0
            Q_ar(i) = 0;
        else
            Q_ar(i) = Q_ar(i) + forces(find(forces(:,1)==i),2);
        end
    elseif isempty(displacements) == 0
        if ismember(i,displacements(:,1))
            Q_ar(i) = NaN;
            arr_unknown_Q = [arr_unknown_Q,i];
        end
    end
end

alpha = arr_unknown_Q;
beta = arr_unknown_q;

qa = q_ar(alpha);
Qb = Q_ar(beta);
qb = K(beta,beta)\(Qb - K(beta,alpha)*qa);
Qa = K(alpha, alpha)*qa + K(alpha,beta)*qb;

Q_sol(alpha) = Qa;
Q_sol(beta) = Qb;
q_sol(alpha) = qa;
q_sol(beta) = qb;

zero_arr = zeros(10,1);
vx_zarray = sym(zero_arr);
vx_yarray = sym(zero_arr);
axial_disp_array = sym(zero_arr);

sheary_nodes = [Q_sol(2),Q_sol(7),Q_sol(12),Q_sol(17),Q_sol(22),Q_sol(27),Q_sol(32),Q_sol(37),Q_sol(42),Q_sol(47),Q_sol(52)];
shearz_nodes = [Q_sol(4),Q_sol(9),Q_sol(14),Q_sol(19),Q_sol(24),Q_sol(29),Q_sol(34),Q_sol(39),Q_sol(44),Q_sol(49),Q_sol(54)];
momentz_nodes = [Q_sol(3),Q_sol(8),Q_sol(13),Q_sol(18),Q_sol(23),Q_sol(28),Q_sol(33),Q_sol(38),Q_sol(43),Q_sol(48),Q_sol(53)];
momenty_nodes = [Q_sol(5),Q_sol(10),Q_sol(15),Q_sol(20),Q_sol(25),Q_sol(30),Q_sol(35),Q_sol(40),Q_sol(45),Q_sol(50),Q_sol(55)];

save('Shear_Moments.mat','nodeloc','sheary_nodes','shearz_nodes','momentz_nodes','momenty_nodes')

%METRIC OR IMPERIAL?
metric = 0; %metric = 1, plot metric. metric = 0, plot imperial

if metric == 1
    figure(1)
    title("Axial Displacement vs Streamwise Distance")
    xlabel("Streamwise Distance Along the Vehicle (m)")
    ylabel("Axial Displacement (m)")
    
    figure(2)
    title("Bending Displacement (z-direction) vs Streamwise Distance")
    xlabel("Streamwise Distance Along the Vehicle (m)")
    ylabel("Bending Displacement (m), z-direction")
    
    figure(3)
    title("Bending Displacement (y-direction) vs Streamwise Distance")
    xlabel("Streamwise Distance Along the Vehicle (m)")
    ylabel("Bending Displacement (m), y-direction")
    
    figure(4)
    title("Axial Stress vs Streamwise Distance")
    xlabel("Streamwise Distance Along the Vehicle (m)")
    ylabel("Axial Stress (Pa)")
    
    figure(5)
    title("Bending Rotation (Around y-axis) vs Streamwise Distance")
    xlabel("Streamwise Distance Along the Vehicle (m)")
    ylabel("Bending Rotation (rad), y-axis")
    
    figure(6)
    title("Bending Moment (Around y-axis)")
    xlabel("Streamwise Distance Along the Vehicle (m)")
    ylabel("Bending Moment (Nm), y-axis")
    
    figure(7)
    title("Bending Shear (z-direction)")
    xlabel("Streamwise Distance Along the Vehicle (m)")
    ylabel("Bending Shear (N) , z-direction")
    
    figure(8)
    title("Bending Rotation (Around z-axis) vs Streamwise Distance")
    xlabel("Streamwise Distance Along the Vehicle (m)")
    ylabel("Bending Rotation (rad), z-axis")
    
    figure(9)
    title("Bending Moment (Around z-axis)")
    xlabel("Streamwise Distance Along the Vehicle (m)")
    ylabel("Bending Moment (Nm), z-axis")
    
    figure(10)
    title("Bending Shear (y-direction)")
    xlabel("Streamwise Distance Along the Vehicle (m)")
    ylabel("Bending Shear (N) , y-direction")

    for i = 1:10
        Nt_local = subs(Nt_x, [E,A,L,I,x], [Eels(i),Aels(i),Lels(i),Iels(i),x-nodeloc(i)]);
        Na_local = subs(Na, [E,A,L,I,x], [Eels(i),Aels(i),Lels(i),Iels(i),x-nodeloc(i)]);
        vx_zarray(i) = q_sol(5*i-3)*Nt_local(1) + q_sol(5*i-2)*Nt_local(2) + q_sol(5*i+2)*Nt_local(3) + q_sol(5*i+3)*Nt_local(4);
        vx_yarray(i) = q_sol(5*i-1)*Nt_local(1) + q_sol(5*i)*Nt_local(2) + q_sol(5*i+4)*Nt_local(3) + q_sol(5*i+5)*Nt_local(4);
        axial_disp_array(i) = q_sol(5*i-4)*Na_local(1) + q_sol(5*i+1)*Na_local(2);
    
        axial_stress(i)= -Eels(i)*diff(axial_disp_array(i),x);
    
        vx_zp(i) = diff(vx_zarray(i),x);
        vx_zpp(i)= diff(vx_zp(i),x);
        vx_zppp(i) = diff(vx_zpp(i),x);
    
        vx_yp(i) = diff(vx_yarray(i),x);
        vx_ypp(i) = diff(vx_yp(i),x);
        vx_yppp(i) = diff(vx_ypp(i),x);
    
        figure(1)
        hold on
        fplot(axial_disp_array(i), [nodeloc(i),nodeloc(i)+Lels(i)])
    
        figure(2)
        hold on
        fplot(vx_zarray(i), [nodeloc(i),nodeloc(i)+Lels(i)])
    
    
        figure(3)
        hold on
        fplot(vx_yarray(i), [nodeloc(i),nodeloc(i)+Lels(i)])
    
        figure(4)
        hold on
        fplot(axial_stress(i),[nodeloc(i),nodeloc(i)+Lels(i)])
    
        figure(5)
        hold on
        fplot(vx_zp(i), [nodeloc(i),nodeloc(i)+Lels(i)])
    
        figure(6)
        hold on
        fplot(-Eels(i)*Iels(i)*vx_zpp(i), [nodeloc(i),nodeloc(i)+Lels(i)])
    
    
        figure(7)
        hold on
        fplot(-Eels(i)*Iels(i)*vx_zppp(i), [nodeloc(i),nodeloc(i)+Lels(i)])
    
        figure(8)
        hold on
        fplot(vx_yp(i),[nodeloc(i),nodeloc(i)+Lels(i)])
    
    
        figure(9)
        hold on
        fplot(-Eels(i)*Iels(i)*vx_ypp(i), [nodeloc(i),nodeloc(i)+Lels(i)])
    
        figure(10)
        hold on
        fplot(-Eels(i)*Iels(i)*vx_yppp(i),[nodeloc(i),nodeloc(i)+Lels(i)])
    
    
    end
elseif metric == 0
    figure(1)
    title("Axial Displacement vs Streamwise Distance")
    xlabel("Streamwise Distance Along the Vehicle (ft)")
    ylabel("Axial Displacement (in)")
    
    figure(2)
    title("Bending Displacement (z-direction) vs Streamwise Distance")
    xlabel("Streamwise Distance Along the Vehicle (ft)")
    ylabel("Bending Displacement (in), z-direction")
    
    figure(3)
    title("Bending Displacement (y-direction) vs Streamwise Distance")
    xlabel("Streamwise Distance Along the Vehicle (ft)")
    ylabel("Bending Displacement (in), y-direction")
    
    figure(4)
    title("Axial Stress vs Streamwise Distance")
    xlabel("Streamwise Distance Along the Vehicle (ft)")
    ylabel("Axial Stress (psi)")
    
    figure(5)
    title("Bending Rotation (Around y-axis) vs Streamwise Distance")
    xlabel("Streamwise Distance Along the Vehicle (ft)")
    ylabel("Bending Rotation (rad), y-axis")
    
    figure(6)
    title("Bending Moment (Around y-axis)")
    xlabel("Streamwise Distance Along the Vehicle (ft)")
    ylabel("Bending Moment (lb*ft), y-axis")
    
    figure(7)
    title("Bending Shear (z-direction)")
    xlabel("Streamwise Distance Along the Vehicle (ft)")
    ylabel("Bending Shear (lb) , z-direction")
    
    figure(8)
    title("Bending Rotation (Around z-axis) vs Streamwise Distance")
    xlabel("Streamwise Distance Along the Vehicle (ft)")
    ylabel("Bending Rotation (rad), z-axis")
    
    figure(9)
    title("Bending Moment (Around z-axis)")
    xlabel("Streamwise Distance Along the Vehicle (ft)")
    ylabel("Bending Moment (lb*ft), z-axis")
    
    figure(10)
    title("Bending Shear (y-direction)")
    xlabel("Streamwise Distance Along the Vehicle (ft)")
    ylabel("Bending Shear (lb) , y-direction")

    for i = 1:10
        Nt_local = subs(Nt_x, [E,A,L,I,x], [Eels(i),Aels(i),Lels(i),Iels(i),x-nodeloc(i)]);
        Na_local = subs(Na, [E,A,L,I,x], [Eels(i),Aels(i),Lels(i),Iels(i),x-nodeloc(i)]);
        vx_zarray(i) = q_sol(5*i-3)*Nt_local(1) + q_sol(5*i-2)*Nt_local(2) + q_sol(5*i+2)*Nt_local(3) + q_sol(5*i+3)*Nt_local(4);
        vx_yarray(i) = q_sol(5*i-1)*Nt_local(1) + q_sol(5*i)*Nt_local(2) + q_sol(5*i+4)*Nt_local(3) + q_sol(5*i+5)*Nt_local(4);
        axial_disp_array(i) = q_sol(5*i-4)*Na_local(1) + q_sol(5*i+1)*Na_local(2);
    
        axial_stress(i)= -Eels(i)*diff(axial_disp_array(i),x);
    
        vx_zp(i) = diff(vx_zarray(i),x);
        vx_zpp(i)= diff(vx_zp(i),x);
        vx_zppp(i) = diff(vx_zpp(i),x);
    
        vx_yp(i) = diff(vx_yarray(i),x);
        vx_ypp(i) = diff(vx_yp(i),x);
        vx_yppp(i) = diff(vx_ypp(i),x);

        xlocs_m = linspace(nodeloc(i),(nodeloc(i)+Lels(i)));
        xlocs_ft = 3.28084*xlocs_m;

        figure(1)
        hold on
        axial_disp_in = 39.3701*subs(axial_disp_array(i),x,xlocs_m);
        plot(xlocs_ft,axial_disp_in);
    
        figure(2)
        hold on
        vz_disp_in = 39.3701*subs(vx_zarray(i),x,xlocs_m);
        plot(xlocs_ft,vz_disp_in);
    
        figure(3)
        hold on
        vy_disp_in = 39.3701*subs(vx_yarray(i),x,xlocs_m);
        plot(xlocs_ft,vy_disp_in);
    
        figure(4)
        hold on
        psi_axial = 0.000145038*subs(axial_stress(i),x,xlocs_m);
        plot(xlocs_ft, psi_axial)
    
        figure(5)
        hold on
        vxzp_rad = subs(vx_zp(i),x,xlocs_m);
        plot(xlocs_ft, vxzp_rad);
    
        figure(6)
        hold on
        My_lbft = (1/1.356)*subs(-Eels(i)*Iels(i)*vx_zpp(i),x,xlocs_m);
        plot(xlocs_ft,My_lbft)
    
        figure(7)
        hold on
        Sz_lb = 0.2248*subs(-Eels(i)*Iels(i)*vx_zppp(i),x,xlocs_m);
        plot(xlocs_ft,Sz_lb);
    
        figure(8)
        hold on
        vxyp_rad = subs(vx_yp(i),x,xlocs_m);
        plot(xlocs_ft, vxyp_rad);
    
        figure(9)
        hold on
        Mz_lbft = (1/1.356)*subs(-Eels(i)*Iels(i)*vx_ypp(i),x,xlocs_m);
        plot(xlocs_ft,Mz_lbft)
    
        figure(10)
        hold on
        Sy_lb = 0.2248*subs(-Eels(i)*Iels(i)*vx_yppp(i),x,xlocs_m);
        plot(xlocs_ft,Sy_lb);
    end
end





function [force] = shearmoment(px, N,Lel,elnode1loc)
    syms x L
    N = subs(N, L, Lel);
    px = subs(px,x, x-elnode1loc);
    F1 = int(px*N(1),x,0,Lel);
    F2 = int(px*N(2),x,0,Lel);
    F3 = int(px*N(3),x,0,Lel);
    F4 = int(px*N(4),x,0,Lel);
    force = [F1,F2,F3,F4];
end
