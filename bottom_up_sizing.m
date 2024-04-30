%Rocket Sizing
%
%variable convension
%   section.part.dims = [length (m), mass (kg), Izz (m^4), apx thickness, CG (relative to forward edge of component), CG (relative to nose of rocket)]
%   section.part.mat = [rho, E, sigY, sigU]
%
%Output:
%   vehicle_sizing.property
%   properties:
%       mass = [first_stage, second_stage,total]
%       lengths = [first_stage,second_stage,total]
%       first_stage = all first stage components and properties
%       second_stage = all second stage components and properties
%       vehicle_size.CS = [front, side]
%
%functions:
%   material(input)
%   dome(aspect_ratio, diameter)
%
%
function [vehicle_sizing] = bottom_up_sizing(solid_prop_mass, liquid_prop_mass, Din)
    
    kg2lb = 2.20462;
    m2ft = 3.28084;

    solid_prop_mass = solid_prop_mass*1.01;
    liquid_prop_mass = liquid_prop_mass*1.027;
    
    %%Props Dimensional Calcs
    %solid propellant & liner
    solid_prop_density = 1763.66;
    
    %solid_prop_density = 1909.8;
    grain_cutout_percent = 10;
    solid_volume = (1+0.01*grain_cutout_percent)*(solid_prop_mass/solid_prop_density);
    [solid_dSA, solid_dH, solid_dV] = dome(sqrt(2), Din);
    solid_cyl_len = (solid_volume-solid_dV)/(pi*(Din/2)^2);
    grain_cutout_radius = sqrt((solid_volume - solid_prop_mass/solid_prop_density)/(pi*(solid_cyl_len + solid_dH)));
    Izz_solid_prop = (solid_prop_mass/2)*((Din/2)^2+grain_cutout_radius^2);
    first_stage.solidprop.dims = [solid_cyl_len+solid_dH, solid_prop_mass, Izz_solid_prop,0, solid_dH+solid_cyl_len/2];

    vehicle_sizing.dH.sol = solid_dH;

    %first_stage_insulation
    insulation_thickness = 0.0127;
    first_stage.insl.mat = material(7); %insulation;
    insl_dome_mass = solid_dSA*insulation_thickness*first_stage.insl.mat(1);
    first_stage.insl.dims = [solid_dH,insl_dome_mass, 0,insulation_thickness,solid_dH/2];

    apx_sol_t = 0.0042; %approximate thickness for solid casing
    [solid2_dSA, solid2_dH, solid2_dV] = dome(sqrt(2), Din+2*insulation_thickness);
    first_stage.solidcasing.mat = material(6); %carbon fiber #2
    sol_casing_cyl_mass = solid_cyl_len*pi*(Din)*apx_sol_t*first_stage.solidcasing.mat(1);
    sol_casing_dome_mass = 2*solid2_dSA*apx_sol_t*first_stage.solidcasing.mat(1);
    sol_casing_mass = sol_casing_cyl_mass+sol_casing_dome_mass;
    first_stage.solidcasing.dims = [solid_cyl_len+2*solid2_dH,sol_casing_mass,(pi/64)*((Din+2*apx_sol_t)^4-(Din)^4), apx_sol_t,solid2_dH+solid_cyl_len/2];
    %solid2_dSA*2 + solid_cyl_len*pi*(Din + 0.012)

    %liquid propellant
    O2F = 1.325;
    %O2F = 1.275;
    oxidizer_mass = (O2F)/(O2F+1)*liquid_prop_mass;
    fuel_mass = 1/(O2F+1)*liquid_prop_mass;
    ox_density = 1440; %N2O4
    fuel_density = 1004.5; %Hydrazine
    ullage_oxpercent =  4;
    ullage_fuelpercent =  4;
    ox_volume = oxidizer_mass/ox_density;
    fuel_volume = fuel_mass/fuel_density;
    %assuming ox on top, and common dome goes concave-up
    [liquid_dSA, liquid_dH, liquid_dV] = dome(sqrt(2), Din);
    if fuel_volume*(1+0.01*ullage_oxpercent) < 2*liquid_dV
        fuel_tank_vol = 2*liquid_dV;
        true_fuel_ullage = (2*liquid_dV/fuel_volume - 1);
        fuel_cyl_len = 0;
    else
        fuel_tank_vol = (1 + 0.01*ullage_fuelpercent)*fuel_volume;
        fuel_cyl_len = (fuel_tank_vol-2*liquid_dV)/(pi*(Din/2)^2);
    end
    ox_tank_vol = (1 + 0.01*ullage_oxpercent)*ox_volume;
    ox_cyl_len = ox_tank_vol/(pi*(Din/2)^2);
    second_stage.ox.dims = [0,oxidizer_mass,0,0,liquid_dH];
    second_stage.fuel.dims = [0,fuel_mass,0,0,liquid_dH + fuel_cyl_len/2];

    %liquid casings
    apx_liq_t = 0.00094;
    second_stage.liqcasing.mat = material(6); %carbon fiber
    liq_tank_len = 2*liquid_dH+ox_cyl_len+fuel_cyl_len+2*apx_liq_t;
    liq_dome_mass = liquid_dSA*apx_liq_t*second_stage.liqcasing.mat(1);
    liq_oxcyl_mass = Din*pi*ox_cyl_len*apx_liq_t*second_stage.liqcasing.mat(1);
    liq_fuelcyl_mass = Din*pi*fuel_cyl_len*apx_liq_t*second_stage.liqcasing.mat(1);
    liq_tank_mass = 3*liq_dome_mass+liq_fuelcyl_mass+liq_oxcyl_mass;
    casing_CoM = ((liq_dome_mass*liquid_dH/2) + liq_dome_mass*(liquid_dH+ox_cyl_len+liquid_dH/2) + liq_oxcyl_mass*(liquid_dH+ox_cyl_len/2) + liq_fuelcyl_mass*(fuel_cyl_len/2 + 2*liquid_dH + ox_cyl_len) + liq_dome_mass*(liquid_dH/2 + fuel_cyl_len + 2*liquid_dH + ox_cyl_len))/liq_tank_mass;
    second_stage.liqcasing.dims = [liq_tank_len, liq_tank_mass, (pi/64)*((Din+2*apx_liq_t)^4-Din^4), apx_liq_t, casing_CoM];

    %aft-skirt and engine
    first_stage.engine.dims = [0.2932176, 7,0,0,0.1];
    apx_aftskirt_t = 0.001;
    aftskirt_len = solid_dH + first_stage.engine.dims(1)/4;
    first_stage.aftskirt.mat = material(1); %aluminum
    apx_fin_mass = 5;
    aftskirt_mass = aftskirt_len*pi*Din*apx_aftskirt_t*first_stage.aftskirt.mat(1) + apx_fin_mass;
    first_stage.aftskirt.dims = [aftskirt_len, aftskirt_mass, (pi/64)*((Din+2*apx_aftskirt_t)^4-Din^4), apx_aftskirt_t, solid_dH + first_stage.engine.dims(1)/8];
    
    first_stage.avgnc.dims = [0.1,1.81437 + 1.67 + 3.62874 + 12.247,0,0,0]; %approximate avionics and gnc mass

    insl_aftskirt_mass = (Din*pi*aftskirt_len)*insulation_thickness*first_stage.insl.mat(1);
    %first_stage.insl.dims = first_stage.insl.dims + [aftskirt_len-solid_dH,insl_aftskirt_mass, 0,0,aftskirt_len/4];

    
    %interstage + wiring runners + liquid engine
    second_stage.engine.dims = [0.5297+0.15, 8,0,0,0.2];
    
    apx_interstage_t = 0.001;
    interstage_len = (liquid_dH + second_stage.engine.dims(1) + solid_dH);
    first_stage.interstage.mat = material(6);
    interstage_mass = interstage_len*pi*Din*apx_interstage_t*first_stage.interstage.mat(1);
    first_stage.interstage.dims = [interstage_len, interstage_mass, (pi/64)*((Din+2*apx_interstage_t)^4-Din^4), apx_interstage_t, interstage_len/2];
    
    second_stage.wiringgnc.dims = [interstage_len/2 + second_stage.liqcasing.dims(1), 1.67 + 2.49476 + 0.453592, 0, 0, 0.5*(interstage_len/2 + second_stage.liqcasing.dims(1))];%FIND IZZ - use cylinders + parallel axis thm
    
    vehicle_sizing.dH.liq = liquid_dH;

    %payload fairing + PAF + payload - finalize when we know final design of PLF
    apx_PLF_t = 0.001;
    first_stage.PLF.dims = [0.5, 5, (pi/64)*((Din+2*apx_PLF_t)^4-Din^4),apx_PLF_t,0.35];
    first_stage.PLF.mat = material(1);%aluminum
    second_stage.PAF.dims = [liquid_dH,3+1,0,0,0.01]; %use plate Izz? 2 kg for MLB, 1.68 for flight computer
    second_stage.payload.dims = [0.3048, 45.3592, 0,0,0.3048/2 ]; %use cube Izz?
    
    %%Full Rocket Mass & CG Calcs
    first_stage_mass = first_stage.solidprop.dims(2) + first_stage.solidcasing.dims(2) + first_stage.engine.dims(2) + first_stage.aftskirt.dims(2) + first_stage.avgnc.dims(2) + first_stage.PLF.dims(2)+ first_stage.interstage.dims(2)*0.8 + first_stage.insl.dims(2);
    second_stage_mass = second_stage.ox.dims(2) + second_stage.fuel.dims(2) + second_stage.liqcasing.dims(2) + second_stage.engine.dims(2) + second_stage.wiringgnc.dims(2) + second_stage.PAF.dims(2) + second_stage.payload.dims(2)+ first_stage.interstage.dims(2)*0.2;
    total_mass = first_stage_mass + second_stage_mass;
    
    first_stage_length = first_stage.engine.dims(1)+first_stage.solidcasing.dims(1);
    second_stage_length = (interstage_len - liquid_dH - solid_dH) + second_stage.liqcasing.dims(1) + first_stage.PLF.dims(1);
    total_length = first_stage_length + second_stage_length;
    
    first_stage_structural_ratio = (first_stage_mass - first_stage.solidprop.dims(2)/1.01)/first_stage_mass;
    second_stage_structural_ratio = (second_stage_mass - (second_stage.ox.dims(2)+second_stage.fuel.dims(2))/1.027 - second_stage.payload.dims(2))/(second_stage_mass - second_stage.payload.dims(2));
    
    first_stage.PLF.dims(6)     =   first_stage.PLF.dims(5);
    second_stage.payload.dims(6)=  second_stage.payload.dims(5)+first_stage.PLF.dims(1) - 1/m2ft;
    second_stage.PAF.dims(6)    =   first_stage.PLF.dims(1) +second_stage.PAF.dims(5);
    second_stage.liqcasing.dims(6)= first_stage.PLF.dims(1)+second_stage.liqcasing.dims(5);
    second_stage.ox.dims(6) =       second_stage.ox.dims(5) + first_stage.PLF.dims(1);
    second_stage.fuel.dims(6) =     second_stage.fuel.dims(5) + second_stage.ox.dims(1) + first_stage.PLF.dims(1);
    first_stage.interstage.dims(6)= first_stage.interstage.dims(5)+second_stage.liqcasing.dims(1) + first_stage.PLF.dims(1)-liquid_dH;
    second_stage.wiringgnc.dims(6)=    second_stage.wiringgnc.dims(5) + first_stage.PLF.dims(1);
    second_stage.engine.dims(6)=    second_stage.engine.dims(5)+second_stage.liqcasing.dims(1) + first_stage.PLF.dims(1);
    first_stage.solidcasing.dims(6)= second_stage_length + first_stage.solidcasing.dims(5);
    first_stage.insl.dims(6)= second_stage_length + first_stage.insl.dims(5);
    first_stage.solidprop.dims(6) = second_stage_length + first_stage.solidprop.dims(5);
    first_stage.aftskirt.dims(6)=   second_stage_length + first_stage.solidcasing.dims(1) - solid_dH + first_stage.aftskirt.dims(5);
    first_stage.engine.dims(6)=     second_stage_length + first_stage.solidcasing.dims(1) + first_stage.engine.dims(5);
    first_stage.avgnc.dims(6)=      second_stage_length + first_stage.solidcasing.dims(1);

    vehicle_CG = (first_stage.PLF.dims(2)*first_stage.PLF.dims(5) +...
        second_stage.payload.dims(2)*second_stage.payload.dims(6)+ ...
        second_stage.PAF.dims(2)*second_stage.PAF.dims(6) +...
        second_stage.liqcasing.dims(2)*second_stage.liqcasing.dims(6) +...
        second_stage.ox.dims(2)*second_stage.ox.dims(6)+...
        second_stage.fuel.dims(2)*second_stage.fuel.dims(6)+...
        first_stage.interstage.dims(2)*first_stage.interstage.dims(6) +...
        second_stage.wiringgnc.dims(2)*second_stage.wiringgnc.dims(6)+...
        second_stage.engine.dims(2)*second_stage.engine.dims(6)+...
        first_stage.solidcasing.dims(2)*first_stage.solidcasing.dims(6)+...
        first_stage.insl.dims(2)*first_stage.insl.dims(6)+...
        first_stage.solidprop.dims(2)*first_stage.solidprop.dims(6)+...
        first_stage.aftskirt.dims(2)*first_stage.aftskirt.dims(6) + ...
        first_stage.engine.dims(2)*first_stage.engine.dims(6) +...
        first_stage.avgnc.dims(2)*first_stage.avgnc.dims(6))/total_mass;
    
    disp(strcat("(metric) Total Vehicle Length: ", num2str(total_length),"m, Vehicle Total Mass: ", num2str(total_mass),"kg"))
    disp(strcat("First Stage - m0: ",num2str(first_stage_mass), "kg, structural ratio: ",num2str(first_stage_structural_ratio*100),"%", ", First Stage Vehicle Length: ", num2str(first_stage_length)))
    disp(strcat("Second Stage - m0: ",num2str(second_stage_mass), "kg, structural ratio: ",num2str(second_stage_structural_ratio*100),"%", ", Second Stage Vehicle Length: ", num2str(second_stage_length)))
    disp(strcat("First Stage Structural Mass: ",num2str(first_stage_structural_ratio*first_stage_mass), "kg, Second Stage Structural Mass: ",num2str(second_stage_structural_ratio*(second_stage_mass-second_stage.payload.dims(2))),"kg"))
    fprintf('\n') 
    disp(strcat("(imperial) Total Vehicle Length: ", num2str(total_length*m2ft),"ft, Vehicle Total Mass: ", num2str(total_mass*kg2lb),"lb"))
    disp(strcat("First Stage - m0: ",num2str(first_stage_mass*kg2lb), "lb, structural ratio: ",num2str(first_stage_structural_ratio*100),"%", ", First Stage Vehicle Length: ", num2str(first_stage_length*m2ft)))
    disp(strcat("Second Stage - m0: ",num2str(second_stage_mass*kg2lb), "lb, structural ratio: ",num2str(second_stage_structural_ratio*100),"%", ", Second Stage Vehicle Length: ", num2str(second_stage_length*m2ft)))
        disp(strcat("First Stage Structural Mass: ",num2str(first_stage_structural_ratio*first_stage_mass*kg2lb), "lb, Second Stage Structural Mass: ",num2str(second_stage_structural_ratio*(second_stage_mass-second_stage.payload.dims(2))*kg2lb),"lb"))
    
    cross_sec_front = pi*(Din/2)^2;
    cross_sec_side = Din*(first_stage_length + second_stage_length);

    vehicle_sizing.first_stage = first_stage;
    vehicle_sizing.second_stage = second_stage;
    vehicle_sizing.CS = [cross_sec_front, cross_sec_side];
    vehicle_sizing.mass = [first_stage_mass, second_stage_mass, total_mass];
    vehicle_sizing.CG = vehicle_CG;
    vehicle_sizing.lengths = [first_stage_length,second_stage_length, total_length];

    save("Vehicle_Sizing.mat",'vehicle_sizing')

    %%Functions
    
    %   material(input)
    %       input:
    %       1 - Aluminum
    %       2 - Stainless Steel
    %       3 - Inconel
    %       4 - Titanium
    %       5 - Carbon Fiber
    
    function [props] = material(input)
        if input == 1 %Aluminum
            rho = 2700;
            E = 7.31E+10;
            sigY = 5.90E+08;
            sigU = 6.41E+08;
        elseif input == 2 %Stainless Steel
            rho = 8000;
            E = 1.93*10^11;
            sigY = 2.05E+08;
            sigU = 5.15E+08;
        elseif input == 3 %Inconel
            rho = 8190;
            E = 2E11;
            sigY = 1.1E+09;
            sigU = 1.38E+09;
        elseif input == 4 %Titanium
            rho = 4430;
            E = 1.14E+11;
            sigY = 8.8E+08;
            sigU = 9.5E+08;
        elseif input == 5 %Carbon Fiber
            rho = 1760;
            E = 2.30E+11;
            sigY = 3530000000;
            sigU = 3530000000;
        elseif input == 6 %Carbon Fiber - Prepreg
            rho = 1486.337198;
            E = 2.895798e+11;
            sigY = 5.688e+9;
            sigU = 5.688e+9;
        elseif input == 7 %Liner
            rho = 1219.005062758356;
            E = 0;
            sigY = 0;
            sigU = 0;
        end
     props = [rho,E,sigY,sigU];
    end


    function [surface_area, height, volume] = dome(AR,D)
        R = D/2;
        height = R/AR;
        p = 1.6075;
        surface_area = 4*pi*((R^p*R^p + 2*R^p*height^p)/3).^(1/p);
        volume = (4*pi*height*R^2/3)/2;
    end
end


