%**************************************************************************
% Filename: FFR_Sim_Macro.m
% Group Name: TW-E
% Date: 04/29/2020
% Description: Script used to simulate Fractional Frequency Reuse (FFR)
% medthods for HetNet applications. FFR algorithm is modeled to
% analyze performance tradeoffs. 
% This file calculates the throughput of the femtocells in the network.
%
% FFR Algorithms:
%                FFR-3SL (Proposed Paper) 
%
%**************************************************************************
close all;

%----------------------------------------------------
% Constants
%----------------------------------------------------
noise_PSD           = -174;             % Noise Power Spectral Density (dBm/Hz)
BER                 = 10^-4;            % Target Bit Error Rate(BER)
alpha               = -1.5/log(BER);	% Constant for target BER
delta_f             = 15e3;             % Subcarrier spacing (Hz)
transmitPower_macro	= 20;               % Macrocell transmit power in Watts
transmitPower_femto = 20e-3;            % Femtocell Base Station Transmit Power in Watts
CH_BW               = 20e6;             % Channel Bandwidth (Hz)
max_subcarriers     = CH_BW/delta_f;    % Maximum number of subcarriers

n_macro_total       = 7;        % Number of macrocells in HetNet
n_femto_per_macro   = 30;       % Number of femtocells per macrocell
n_femto_total       = 210;      % Total femtocells (number_per_macro * 7)
r_macro             = 500;      % Radius of macrocell hexagon
n_mue_per_macro     = 150;      % Number of UEs per macrocell
n_mue_total         = 1050;     % Total UEs in the network
number_of_runs      = 100;      % Number of times to run the Monte Carlo sim

% Store all possible increments of the femtocell count to loop through
femtocell_array     = 1:n_femto_total;

% Initialize the throughput array
throughput_macro_array = zeros(1,n_femto_total);

%----------------------------------------------------
% Set the coordinates of the macrocell towers and populate with femotcells
%----------------------------------------------------
% Capture the center of each hexagon (location of each macrocell tower)
A_center_X = 1500;
A_center_Y = 1500;
B_center_X = 1500;
B_center_Y = 2500;
C_center_X = 2250;
C_center_Y = 2000;
D_center_X = 2250;
D_center_Y = 1000;
E_center_X = 1500;
E_center_Y = 500;
F_center_X = 750;
F_center_Y = 1000;
G_center_X = 750;
G_center_Y = 2000;

macro_X_coords = [A_center_X B_center_X C_center_X D_center_X E_center_X F_center_X G_center_X];
macro_Y_coords = [A_center_Y B_center_Y C_center_Y D_center_Y E_center_Y F_center_Y G_center_Y];

% Define the vertexes of the hexagon for angles 0, 60, 120, 180, 240 and
% 300 with the origin offset by the center X,Y coordinates
Av_x = (r_macro * cos((0:6)*pi/3)) + A_center_X;
Av_y = (r_macro * sin((0:6)*pi/3)) + A_center_Y;
Bv_x = (r_macro * cos((0:6)*pi/3)) + B_center_X;
Bv_y = (r_macro * sin((0:6)*pi/3)) + B_center_Y;
Cv_x = (r_macro * cos((0:6)*pi/3)) + C_center_X;
Cv_y = (r_macro * sin((0:6)*pi/3)) + C_center_Y;
Dv_x = (r_macro * cos((0:6)*pi/3)) + D_center_X;
Dv_y = (r_macro * sin((0:6)*pi/3)) + D_center_Y;
Ev_x = (r_macro * cos((0:6)*pi/3)) + E_center_X;
Ev_y = (r_macro * sin((0:6)*pi/3)) + E_center_Y;
Fv_x = (r_macro * cos((0:6)*pi/3)) + F_center_X;
Fv_y = (r_macro * sin((0:6)*pi/3)) + F_center_Y;
Gv_x = (r_macro * cos((0:6)*pi/3)) + G_center_X;
Gv_y = (r_macro * sin((0:6)*pi/3)) + G_center_Y;

% Generate many points in a square and choose N points that fall within the hexagon
% Generate 3*n_femto_per_macro random points with square that is 2R by 2R
Ac_x = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + A_center_X;
Ac_y = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + A_center_Y;
Bc_x = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + B_center_X;
Bc_y = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + B_center_Y;
Cc_x = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + C_center_X;
Cc_y = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + C_center_Y;
Dc_x = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + D_center_X;
Dc_y = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + D_center_Y;
Ec_x = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + E_center_X;
Ec_y = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + E_center_Y;
Fc_x = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + F_center_X;
Fc_y = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + F_center_Y;
Gc_x = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + G_center_X;
Gc_y = (r_macro-rand(1, 3*n_femto_per_macro)*2*r_macro) + G_center_Y;

% The MATLAB command inploygon finds points within a polygon region.
A_IN = inpolygon(Ac_x, Ac_y, Av_x, Av_y);
B_IN = inpolygon(Bc_x, Bc_y, Bv_x, Bv_y);
C_IN = inpolygon(Cc_x, Cc_y, Cv_x, Cv_y);
D_IN = inpolygon(Dc_x, Dc_y, Dv_x, Dv_y);
E_IN = inpolygon(Ec_x, Ec_y, Ev_x, Ev_y);
F_IN = inpolygon(Fc_x, Fc_y, Fv_x, Fv_y);
G_IN = inpolygon(Gc_x, Gc_y, Gv_x, Gv_y);

% Drop femtocell nodes outside the hexagon
Ac_x = Ac_x(A_IN);
Ac_y = Ac_y(A_IN);
Bc_x = Bc_x(B_IN);
Bc_y = Bc_y(B_IN);
Cc_x = Cc_x(C_IN);
Cc_y = Cc_y(C_IN);
Dc_x = Dc_x(D_IN);
Dc_y = Dc_y(D_IN);
Ec_x = Ec_x(E_IN);
Ec_y = Ec_y(E_IN);
Fc_x = Fc_x(F_IN);
Fc_y = Fc_y(F_IN);
Gc_x = Gc_x(G_IN);
Gc_y = Gc_y(G_IN);

% Choose only n_femto_per_macro points for femtocells
A_idx = randperm(length(Ac_x));
Ac_x = Ac_x(A_idx(1:n_femto_per_macro));
Ac_y = Ac_y(A_idx(1:n_femto_per_macro));
B_idx = randperm(length(Bc_x));
Bc_x = Bc_x(B_idx(1:n_femto_per_macro));
Bc_y = Bc_y(B_idx(1:n_femto_per_macro));
C_idx = randperm(length(Cc_x));
Cc_x = Cc_x(C_idx(1:n_femto_per_macro));
Cc_y = Cc_y(C_idx(1:n_femto_per_macro));
D_idx = randperm(length(Dc_x));
Dc_x = Dc_x(D_idx(1:n_femto_per_macro));
Dc_y = Dc_y(D_idx(1:n_femto_per_macro));
E_idx = randperm(length(Ec_x));
Ec_x = Ec_x(E_idx(1:n_femto_per_macro));
Ec_y = Ec_y(E_idx(1:n_femto_per_macro));
F_idx = randperm(length(Fc_x));
Fc_x = Fc_x(F_idx(1:n_femto_per_macro));
Fc_y = Fc_y(F_idx(1:n_femto_per_macro));
G_idx = randperm(length(Gc_x));
Gc_x = Gc_x(G_idx(1:n_femto_per_macro));
Gc_y = Gc_y(G_idx(1:n_femto_per_macro));

% Combine all femto coordinates from every macrocell into one array
femto_X_coords = [Ac_x Bc_x Cc_x Dc_x Ec_x Fc_x Gc_x];
femto_Y_coords = [Ac_y Bc_y Cc_y Dc_y Ec_y Fc_y Gc_y];

% Increment through all of the possible femtocell increments
for total_femto_count = 10:10:n_femto_total
    
    % hold the sum of all FUE throughput - reset at every femtocell
    % increment
    total_throughput = 0;
    
    % Initialize the subcarrier assignment counter to 0
    subcarriers_assigned = 0;
    
    % run Monte Carlo Simulation 100 times at every increment of femtocell
    for i = 1:number_of_runs
    
        %--------- Generate new locations for FUEs ----------------------%
        % Generate many points in a square and choose N points that fall within the hexagon
        % Generate 3*n_fue_per_macro random points with square that is 2R by 2R
        A_mue_x = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + A_center_X;
        A_mue_y = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + A_center_Y;
        B_mue_x = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + B_center_X;
        B_mue_y = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + B_center_Y;
        C_mue_x = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + C_center_X;
        C_mue_y = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + C_center_Y;
        D_mue_x = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + D_center_X;
        D_mue_y = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + D_center_Y;
        E_mue_x = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + E_center_X;
        E_mue_y = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + E_center_Y;
        F_mue_x = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + F_center_X;
        F_mue_y = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + F_center_Y;
        G_mue_x = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + G_center_X;
        G_mue_y = (r_macro-rand(1, 3*n_mue_per_macro)*2*r_macro) + G_center_Y;

        % The MATLAB command inploygon finds points within a polygon region
        A_IN_mue = inpolygon(A_mue_x, A_mue_y, Av_x, Av_y);
        B_IN_mue = inpolygon(B_mue_x, B_mue_y, Bv_x, Bv_y);
        C_IN_mue = inpolygon(C_mue_x, C_mue_y, Cv_x, Cv_y);
        D_IN_mue = inpolygon(D_mue_x, D_mue_y, Dv_x, Dv_y);
        E_IN_mue = inpolygon(E_mue_x, E_mue_y, Ev_x, Ev_y);
        F_IN_mue = inpolygon(F_mue_x, F_mue_y, Fv_x, Fv_y);
        G_IN_mue = inpolygon(G_mue_x, G_mue_y, Gv_x, Gv_y);

        % Drop FUE nodes outside the hexagon
        A_mue_x = A_mue_x(A_IN_mue);
        A_mue_y = A_mue_y(A_IN_mue);
        B_mue_x = B_mue_x(B_IN_mue);
        B_mue_y = B_mue_y(B_IN_mue);
        C_mue_x = C_mue_x(C_IN_mue);
        C_mue_y = C_mue_y(C_IN_mue);
        D_mue_x = D_mue_x(D_IN_mue);
        D_mue_y = D_mue_y(D_IN_mue);
        E_mue_x = E_mue_x(E_IN_mue);
        E_mue_y = E_mue_y(E_IN_mue);
        F_mue_x = F_mue_x(F_IN_mue);
        F_mue_y = F_mue_y(F_IN_mue);
        G_mue_x = G_mue_x(G_IN_mue);
        G_mue_y = G_mue_y(G_IN_mue);

        % Choose only n_fue_per_macro points for FUEs
        A_idx = randperm(length(A_mue_x));
        A_mue_x = A_mue_x(A_idx(1:n_mue_per_macro));
        A_mue_y = A_mue_y(A_idx(1:n_mue_per_macro));
        B_idx = randperm(length(B_mue_x));
        B_mue_x = B_mue_x(B_idx(1:n_mue_per_macro));
        B_mue_y = B_mue_y(B_idx(1:n_mue_per_macro));
        C_idx = randperm(length(C_mue_x));
        C_mue_x = C_mue_x(C_idx(1:n_mue_per_macro));
        C_mue_y = C_mue_y(C_idx(1:n_mue_per_macro));
        D_idx = randperm(length(D_mue_x));
        D_mue_x = D_mue_x(D_idx(1:n_mue_per_macro));
        D_mue_y = D_mue_y(D_idx(1:n_mue_per_macro));
        E_idx = randperm(length(E_mue_x));
        E_mue_x = E_mue_x(E_idx(1:n_mue_per_macro));
        E_mue_y = E_mue_y(E_idx(1:n_mue_per_macro));
        F_idx = randperm(length(F_mue_x));
        F_mue_x = F_mue_x(F_idx(1:n_mue_per_macro));
        F_mue_y = F_mue_y(F_idx(1:n_mue_per_macro));
        G_idx = randperm(length(G_mue_x));
        G_mue_x = G_mue_x(G_idx(1:n_mue_per_macro));
        G_mue_y = G_mue_y(G_idx(1:n_mue_per_macro));
        
        % Combine all MUE coordinates from every macrocell into one array
        mue_X_coords = [A_mue_x B_mue_x C_mue_x D_mue_x E_mue_x F_mue_x G_mue_x];
        mue_Y_coords = [A_mue_y B_mue_y C_mue_y D_mue_y E_mue_y F_mue_y G_mue_y];
        %----------------------------------------------------------------%        
        
        % Loop through all MUEs in the network
        for mue_index = 1:n_mue_total
            %-----------Find the nearest macrocell to assign the MUE to------%
            % Initialize shortest distance as the diameter of the network
            d_shortest_macro = 3000;
            assigned_macrocell_index = 10;
            
            % Loop through all femtocells in this iteration
            for macro_index = 1:n_macro_total
                % create a 2x2 matrix of current MUE coordinates and the
                % coordinates of the current femtocell index
                coord_matrix = [mue_X_coords(mue_index),mue_Y_coords(mue_index);macro_X_coords(macro_index),macro_Y_coords(macro_index)];

                % calculate the distance between those two points
                d_macro = pdist(coord_matrix,'euclidean');

                % if the current distance is less than the saved min distance
                if d_macro <= d_shortest_macro
                    % replace the new shortest distance
                    d_shortest_macro = d_macro;

                    % save the index of the closest macrocell for later use
                    assigned_macrocell_index = macro_index;               
                end
            end
            %----------------------------------------------------------------%            
            
            % Start the calculation of Equation (4)
            % multiply noise spectral density by the subcarrier spacing and convert
            % to Watts
            denominator = 10^((noise_PSD * delta_f)/10);

            % Calculate macrocell interference
            % Summation of neighboring macrocell's power & gain products on subcarrier k
            sigma_Pkm_GkmM = 0; % Initialize to zero
            
            % Calculate interference to MUE from interferer macrocells
            for current_macro_index = 1:n_macro_total
                
                % Don't calculate the interference of the assigned macrocell
                if not(current_macro_index == assigned_macrocell_index)
                    % create a 2x2 matrix of current MUE coordinates and the
                    % coordinates of the interfering macrocell. 
                    coord_matrix = [mue_X_coords(mue_index),mue_Y_coords(mue_index);macro_X_coords(current_macro_index),macro_Y_coords(current_macro_index)];

                    % calculate the distance between those two points
                    d_macro = pdist(coord_matrix,'euclidean');

                    % Outdoor Pathloss - Equation (2) from paper
                    PL_macro = 28.0 + 35*log10(d_macro);

                    % Channel Gain - Equation (3) from paper
                    CG_macro = 10^-(PL_macro/10);

                    % Add up all the interferers
                    sigma_Pkm_GkmM = sigma_Pkm_GkmM + (transmitPower_macro*CG_macro);
                end
            end
            
            % Add the macrocell interferers to the denominator of Equation
            % (4)
            denominator = denominator + sigma_Pkm_GkmM;

            % Femtocell interference
            % Summation of F neighboring femtocell power & gain products on sub-carrier k
            sigma_PkF_GkmF = 0; % Initialize to zero

            % Loop through all femtocells to calculate interference
            for current_f_index = 1:total_femto_count           
                % Create a 2x2 matrix of current FUE coordinates and the
                % Coordinates of the interfering femtocell. 
                coord_matrix = [mue_X_coords(mue_index),mue_Y_coords(mue_index);femto_X_coords(current_f_index),femto_Y_coords(current_f_index)];

                % Calculate the distance between those two points
                d_femto = pdist(coord_matrix,'euclidean');

                % Outdoor Pathloss - Equation (2) from paper
                PL_femto = 28.0 + 35*log10(d_femto);

                % Channel Gain - Equation (3) from paper
                CG_femto = 10^-(PL_femto/10);

                sigma_PkF_GkmF = sigma_PkF_GkmF + (transmitPower_femto*CG_femto);
            end

            % Add the femtocell interferers to the denominator of Equation
            % (4)
            denominator = denominator + sigma_PkF_GkmF;

            %--------Numerator of Equation (4)---------------------------%
            
            % Create a 2x2 matrix of current MUE coordinates and the
            % coordinates of the host femtocell
            coord_matrix = [mue_X_coords(mue_index),mue_Y_coords(mue_index);macro_X_coords(assigned_macrocell_index),macro_Y_coords(assigned_macrocell_index)];

            % Calculate the distance between those two points
            d_fue = pdist(coord_matrix,'euclidean');

            % Outdoor Pathloss - Equation (2) from paper
            PL_fue = 28.0 + 35*log10(d_fue);

            % Channel Gain - Equation (3) from paper
            CG_fue = 10^-(PL_fue/10);

            % Numerator of the SINR equation - Equation (4)
            numerator = transmitPower_femto * CG_fue;

            % Combine values into SINR to complete Equation (4)
            SINR_kf = numerator / denominator;

            % Calculate channel capacity - Equation (6) from paper
            channelCapacity_macro = delta_f * log2(1 + (alpha * SINR_kf));

            % Create subcarrier assignments, but only assign up to the
            % alotted maximum number based on bandwidth and channel spacing
            % Equation (8) from paper
            if subcarriers_assigned < (max_subcarriers - 12)
                throughput_macro = channelCapacity_macro * 12;
                subcarriers_assigned = subcarriers_assigned + 12;
            else
                throughput_macro = 0;
            end
            
            % Summation of the throughput over all MUEs and
            % subcarriers
            total_throughput = total_throughput + throughput_macro;
        end
    end
        
    % Assign the throughput for this macrocell increment to the array
    % Average it over the number of runs
    throughput_macro_array(total_femto_count) = (total_throughput/number_of_runs);
end


%----------------------------------------------------
% Plots
%----------------------------------------------------
figure; 
plot(femtocell_array,throughput_macro_array/1e6, 'x');
xlabel('Number of Femtocells');
ylabel('Throughput (Mbps)');
title('Throughput of the UEs connected with Macrocell');