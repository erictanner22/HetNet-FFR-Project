%**************************************************************************
% Filename: FFR_Simulator.m
% Group Name: GTW-E
% Date: 03/17/2020
% Description: Script used to simulate Fractional Frequency Reuse(FFR)
% medthods for HetNet applications. FFR algorithm is modeled to
% analyze performance tradeoffs. 
%
% FFR Algorithms:
%                FFR-3SL (Proposed Paper) 
%
%**************************************************************************
clear all;
close all;

%----------------------------------------------------
% Common Sim Variables
%----------------------------------------------------
BER     = 10^-4;                % Target Bit Error Rate(BER)
alpha   = -1.5/log(BER);        % Constant for target BER
delta_f  = 15e3;                % Subcarrier spacing (Hz)
CH_BW   = 20e6;                 % Channel Bandwidth (Hz)
Max_Subcarriers = CH_BW/delta_f;% Maximum number of subcarriers
macro_radius = 500;             % Macrocell radius(m)
femto_radius = 30;              % femtocell radius(m)
m_users = 150;                  % Macrocell Users
f_users = 150;                  % Femtocell Users
SCpRB   = 12;                   % Sub-carriers per resource block
Num_RB  = 100;                  % Total number of resource blocks
Num_SC  = Num_RB*SCpRB;         % Number of subcarriers
transmitPower_macro = 20;       % Macrocell transmit power in Watts
MC_TxP_vec  = [10 15 20];       % Macrocell Base Station Transmit Power
transmitPower_femto = 20e-3;	% Femtocell Base Station Transmit Power in Watts
transmitPower_mue   = 1;        % MUE transmit power (ASSUMPTION because the paper doesn't cite one)
Noise_PSD  = -174;              % Noise Power Spectral Density (dBm/Hz)
Num_Mc  = 7;                    % Number of Macrocells
Num_Fc  = 30;                   % Number of Femtocells per macrocell
T_Num_Fc  = 210;                % Total Number of Femtocells
Num_RB  = 100;                  % Number of Resource Blocks
SB_Sect = 4;                    % Subbands per sector
total_subbands = 7;             % Total freq band is divided into 7 subbands A-G
d_vec  = 5;                     % Distance from base station to user in meters
MUE_vec = 1:150;                % vector of macro user equipments to loop through
subcarriers_vec = 1:1333;       % vector of subcarriers. This was calculated by 
                                % dividing the channel bandwidth of 20 MHz by the 
                                % subcarrier spacing of 15 kHz. 
number_of_runs = 100;           % number of times to run the Monte Carlo sim


%**************************************************************************
% Macrocell - SINR, Capacity, and Throughput
%**************************************************************************
% Precompute femtocell distance vector for simulation
% Distance to any interferring femtocell will be 60 - 470m
d_femto_vec = round(rand(1,Num_Fc)*(470) + (2*femto_radius));   

% Precompute Macro User Equipment (MUE) distance vector for simulation
d_MUE_vec  = round(rand(1,m_users)*(macro_radius)+1);

idx=0;
Tm_macro_vec = zeros(1,length(Num_Mc));
% Loop through all macro cells placing 150 UEs and calculating throughput
% for all 7 macrocells.
for Nm=1:Num_Mc
       
    SubCarriers_Assigned = 0;
    Tm_user_vec = zeros(1,m_users);
    for m=1:m_users
        
        % Convert Macrocell power to watts
        MC_TxP = MC_TxP_vec(3);
        MC_TxP_W = 10^(MC_TxP/10);
        
        % Convert Femtocell power to watts
        FC_TxP_W = 10^(transmitPower_femto/10);

        % Summation of M neighboring Macro-cell's Power & Gain products on sub-carrier k
        % Equation 4 - Denomonator middle summation
        sigma_PMp_GMp = 0; % Initialize to zero
        sigma_PF_GF   = 0; % Initialize to zero
        for mc=1:(Num_Mc-1)
            %d = d_vec(m);
            d=866;
            PL_macro = 28.0 + 35*log10(d);
            Gain = 10^-(PL_macro/10);
            sigma_PMp_GMp = sigma_PMp_GMp + (MC_TxP_W*Gain);
        end
        
        % Summation of F neighboring Femto-cell Power & Gain products on sub-carrier k
        for fc=1:(Num_Fc)
            PL_femto = 28.0 + 35*log10(d_femto_vec(fc));
            Gain = 10^-(PL_femto/10);
            sigma_PF_GF = sigma_PF_GF + (FC_TxP_W*Gain);
        end
        
        
        % Calculate channel gain (NOTE: Removing the Xsigma and the |H|
        % Rayleigh Gaussian distribution).
        PL_MUE = 28.0 + 35*log10(d_MUE_vec(m));
        Ch_Gain_W = 10^(-PL_MUE/10);
        
        % SINR equation for a given Macro-cell on sub-carrier k
        % NOTE: this is equation 4 from the paper
        SINRmk = (MC_TxP_W*Ch_Gain_W)/(10^((Noise_PSD*delta_f)/10) + sigma_PMp_GMp + sigma_PF_GF);
        
        % Capacity of macro user m on sub-carrier k
        Cmk = delta_f*log2(1+alpha*SINRmk);
        

        % Calculate Throughput of the Macro-cell across all m users and all
        % subcarriers available to user.
        Tm = 0; % Initialize to zero
        for k=1:SCpRB
            
            % Create subcarrier assignments, but only assign up to the
            % alotted maximum number based on bandwidth and channel spacing
            if SubCarriers_Assigned < Max_Subcarriers
                Beta_km = 1;
                SubCarriers_Assigned = SubCarriers_Assigned + 1;
            else
                Beta_km = 0;
            end
        
            Tm = Tm + Cmk * Beta_km;
        end
        
        % Assign throughput based on macrouser
        Tm_user_vec(m) = Tm;
        
    end
    
    idx = idx+1;
    Tm_macro_vec(idx) = sum(Tm_user_vec);
    
end



%**************************************************************************
% Figures/Plots
%**************************************************************************
figure; 
plot((1:Num_Mc)*30,(Tm_macro_vec/1e6), 'o');
xlabel('Number of Femto-cells');
ylabel('Throughput (Mbps)');
title('Macrocell Throughput');





    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vvv Teresa's Additions vvv

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

%In the code, I will create a hexagon centered at (1500,1500) with radius R. The snipplets can be used in mobile capacity predicts and general systems level simulation of cellular networks.
n_femto = 30;      % Number of femtocells
r_macro = 500;     % Radius of Hexagon
n_mue   = 150;     % number of UEs per macrocell

%Define the vertexes of the hexagon. They for angles 0, 60, 120, 180, 240 and 300 withe origin.
%Vertexes
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

%The method used here is to generate many points in a square and choose N points that fall within the hexagon
%Generate 3*n_femto random points with square that is 2R by 2R
Ac_x = (r_macro-rand(1, 3*n_femto)*2*r_macro) + A_center_X;
Ac_y = (r_macro-rand(1, 3*n_femto)*2*r_macro) + A_center_Y;
Bc_x = (r_macro-rand(1, 3*n_femto)*2*r_macro) + B_center_X;
Bc_y = (r_macro-rand(1, 3*n_femto)*2*r_macro) + B_center_Y;
Cc_x = (r_macro-rand(1, 3*n_femto)*2*r_macro) + C_center_X;
Cc_y = (r_macro-rand(1, 3*n_femto)*2*r_macro) + C_center_Y;
Dc_x = (r_macro-rand(1, 3*n_femto)*2*r_macro) + D_center_X;
Dc_y = (r_macro-rand(1, 3*n_femto)*2*r_macro) + D_center_Y;
Ec_x = (r_macro-rand(1, 3*n_femto)*2*r_macro) + E_center_X;
Ec_y = (r_macro-rand(1, 3*n_femto)*2*r_macro) + E_center_Y;
Fc_x = (r_macro-rand(1, 3*n_femto)*2*r_macro) + F_center_X;
Fc_y = (r_macro-rand(1, 3*n_femto)*2*r_macro) + F_center_Y;
Gc_x = (r_macro-rand(1, 3*n_femto)*2*r_macro) + G_center_X;
Gc_y = (r_macro-rand(1, 3*n_femto)*2*r_macro) + G_center_Y;

%The method used here is to generate many points in a square and choose N points that fall within the hexagon
%Generate 3*n_mue random points with square that is 2R by 2R
A_mue_x = (r_macro-rand(1, 3*n_mue)*2*r_macro) + A_center_X;
A_mue_y = (r_macro-rand(1, 3*n_mue)*2*r_macro) + A_center_Y;
B_mue_x = (r_macro-rand(1, 3*n_mue)*2*r_macro) + B_center_X;
B_mue_y = (r_macro-rand(1, 3*n_mue)*2*r_macro) + B_center_Y;
C_mue_x = (r_macro-rand(1, 3*n_mue)*2*r_macro) + C_center_X;
C_mue_y = (r_macro-rand(1, 3*n_mue)*2*r_macro) + C_center_Y;
D_mue_x = (r_macro-rand(1, 3*n_mue)*2*r_macro) + D_center_X;
D_mue_y = (r_macro-rand(1, 3*n_mue)*2*r_macro) + D_center_Y;
E_mue_x = (r_macro-rand(1, 3*n_mue)*2*r_macro) + E_center_X;
E_mue_y = (r_macro-rand(1, 3*n_mue)*2*r_macro) + E_center_Y;
F_mue_x = (r_macro-rand(1, 3*n_mue)*2*r_macro) + F_center_X;
F_mue_y = (r_macro-rand(1, 3*n_mue)*2*r_macro) + F_center_Y;
G_mue_x = (r_macro-rand(1, 3*n_mue)*2*r_macro) + G_center_X;
G_mue_y = (r_macro-rand(1, 3*n_mue)*2*r_macro) + G_center_Y;

%There is a command in MATLAB inploygon. 
%The command finds points within a polygon region.
%get the femtocell points within the polygon
A_IN = inpolygon(Ac_x, Ac_y, Av_x, Av_y);
B_IN = inpolygon(Bc_x, Bc_y, Bv_x, Bv_y);
C_IN = inpolygon(Cc_x, Cc_y, Cv_x, Cv_y);
D_IN = inpolygon(Dc_x, Dc_y, Dv_x, Dv_y);
E_IN = inpolygon(Ec_x, Ec_y, Ev_x, Ev_y);
F_IN = inpolygon(Fc_x, Fc_y, Fv_x, Fv_y);
G_IN = inpolygon(Gc_x, Gc_y, Gv_x, Gv_y);

%There is a command in MATLAB inploygon. 
%The command finds points within a polygon region.
%get the MUE points within the polygon
A_IN_mue = inpolygon(A_mue_x, A_mue_y, Av_x, Av_y);
B_IN_mue = inpolygon(B_mue_x, B_mue_y, Bv_x, Bv_y);
C_IN_mue = inpolygon(C_mue_x, C_mue_y, Cv_x, Cv_y);
D_IN_mue = inpolygon(D_mue_x, D_mue_y, Dv_x, Dv_y);
E_IN_mue = inpolygon(E_mue_x, E_mue_y, Ev_x, Ev_y);
F_IN_mue = inpolygon(F_mue_x, F_mue_y, Fv_x, Fv_y);
G_IN_mue = inpolygon(G_mue_x, G_mue_y, Gv_x, Gv_y);

%drop femto nodes outside the hexagon
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

%drop MUE nodes outside the hexagon
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

%choose only N_femto points for femtocells
A_idx = randperm(length(Ac_x));
Ac_x = Ac_x(A_idx(1:n_femto));
Ac_y = Ac_y(A_idx(1:n_femto));
B_idx = randperm(length(Bc_x));
Bc_x = Bc_x(B_idx(1:n_femto));
Bc_y = Bc_y(B_idx(1:n_femto));
C_idx = randperm(length(Cc_x));
Cc_x = Cc_x(C_idx(1:n_femto));
Cc_y = Cc_y(C_idx(1:n_femto));
D_idx = randperm(length(Dc_x));
Dc_x = Dc_x(D_idx(1:n_femto));
Dc_y = Dc_y(D_idx(1:n_femto));
E_idx = randperm(length(Ec_x));
Ec_x = Ec_x(E_idx(1:n_femto));
Ec_y = Ec_y(E_idx(1:n_femto));
F_idx = randperm(length(Fc_x));
Fc_x = Fc_x(F_idx(1:n_femto));
Fc_y = Fc_y(F_idx(1:n_femto));
G_idx = randperm(length(Gc_x));
Gc_x = Gc_x(G_idx(1:n_femto));
Gc_y = Gc_y(G_idx(1:n_femto));

%choose only n_mue points for MUEs
A_idx = randperm(length(A_mue_x));
A_mue_x = A_mue_x(A_idx(1:n_mue));
A_mue_y = A_mue_y(A_idx(1:n_mue));
B_idx = randperm(length(B_mue_x));
B_mue_x = B_mue_x(B_idx(1:n_mue));
B_mue_y = B_mue_y(B_idx(1:n_mue));
C_idx = randperm(length(C_mue_x));
C_mue_x = C_mue_x(C_idx(1:n_mue));
C_mue_y = C_mue_y(C_idx(1:n_mue));
D_idx = randperm(length(D_mue_x));
D_mue_x = D_mue_x(D_idx(1:n_mue));
D_mue_y = D_mue_y(D_idx(1:n_mue));
E_idx = randperm(length(E_mue_x));
E_mue_x = E_mue_x(E_idx(1:n_mue));
E_mue_y = E_mue_y(E_idx(1:n_mue));
F_idx = randperm(length(F_mue_x));
F_mue_x = F_mue_x(F_idx(1:n_mue));
F_mue_y = F_mue_y(F_idx(1:n_mue));
G_idx = randperm(length(G_mue_x));
G_mue_x = G_mue_x(G_idx(1:n_mue));
G_mue_y = G_mue_y(G_idx(1:n_mue));


% Plot all of the femtocell locations
plot(Ac_x, Ac_y, 'r*');
hold on;
plot(Bc_x, Bc_y, 'r*');
hold on;
plot(Cc_x, Cc_y, 'r*');
hold on;
plot(Dc_x, Dc_y, 'r*');
hold on;
plot(Ec_x, Ec_y, 'r*');
hold on;
plot(Fc_x, Fc_y, 'r*');
hold on;
plot(Gc_x, Gc_y, 'r*');
hold on;

% Plot all of the MUE locations
plot(A_mue_x, A_mue_y, 'bo');
hold on;
plot(B_mue_x, B_mue_y, 'bo');
hold on;
plot(C_mue_x, C_mue_y, 'bo');
hold on;
plot(D_mue_x, D_mue_y, 'bo');
hold on;
plot(E_mue_x, E_mue_y, 'bo');
hold on;
plot(F_mue_x, F_mue_y, 'bo');
hold on;
plot(G_mue_x, G_mue_y, 'bo');
hold on;

% Plot the hexagon outlines
plot(Av_x,Av_y);
plot(Bv_x,Bv_y);
plot(Cv_x,Cv_y);
plot(Dv_x,Dv_y);
plot(Ev_x,Ev_y);
plot(Fv_x,Fv_y);
plot(Gv_x,Gv_y);

hold on;

axis square;

hold off

% combine all femto coordinates from every macrocell into one array
femto_X_coords = [Ac_x Bc_x Cc_x Dc_x Ec_x Fc_x Gc_x];
femto_Y_coords = [Ac_y Bc_y Cc_y Dc_y Ec_y Fc_y Gc_y];

% combine all MUE coordinates from every macrocell into one array
mue_X_coords = [A_mue_x B_mue_x C_mue_x D_mue_x E_mue_x F_mue_x G_mue_x];
mue_Y_coords = [A_mue_y B_mue_y C_mue_y D_mue_y E_mue_y F_mue_y G_mue_y];

% 1x30 array for holding all throughput values - intialize all to 0
throughput_macro_array = zeros(1,210);
mue_distances   = zeros(1,1050);
femtocell_array = 1:210;
number_of_runs = 1;

% increment total femtocells for graphing
% Only include femtocells within macrocell A
for total_femto_count = femtocell_array

    % run Monte Carlo Simulation 100 times at every increment of femtocell
    for i = 1:number_of_runs
    
        % Loop through all MUEs
        for mue_index = 1:n_mue

            % hold the sum of all MUE throughput - reset at every femtocell
            % increment
            total_throughput = 0;

            % multiply noise spectral density by the subcarrier spacing and convert
            % to Watts
            denominator = 10^((Noise_PSD * delta_f)/10);

            % calculate macrocell interference
            % Summation of M neighboring Macro-cell's Power & Gain products on sub-carrier k
            sigma_Pkm_GkmM = 0; % Initialize to zero

            % Calculate interference to MUE from interferer macrocells
            for current_m_index = 1:Num_Mc
                % figure out which index to skip (don't count the macrocell the
                % MUE is assigned to
                % calculate the starting and ending indices for the MUE array
                % A(1-150), B(151-300), C(301-450), D(451-600), E(601-750),
                % F(751-900), G(901-1050)
                first_mue_index = (current_m_index - 1)*n_mue + 1;
                last_mue_index  = current_m_index*n_mue;

                % Only do the SINR calc for macrocells that aren't the desired
                % one
                if not(mue_index >= first_mue_index && mue_index <= last_mue_index)
                    % create a 2x2 matrix of current mue coordinates and the
                    % coordinates of the interfering macrocell. 
                    coord_matrix = [mue_X_coords(mue_index),mue_Y_coords(mue_index);macro_X_coords(current_m_index),macro_Y_coords(current_m_index)];

                    % calculate the distance between those two points
                    d_macro = pdist(coord_matrix,'euclidean');

                    % outdoor pathloss - equation (2) from paper
                    PL_macro = 28.0 + 35*log10(d_macro);

                    % equation (3) from paper
                    CG_macro = 10^-(PL_macro/10);

                    % Add up all the interferers
                    sigma_Pkm_GkmM = sigma_Pkm_GkmM + (transmitPower_macro*CG_macro);
                end
            end

            % add the macrocell interferers to the denom 
            denominator = denominator + sigma_Pkm_GkmM;

            % femtocell interference
            % Summation of F neighboring Femto-cell Power & Gain products on sub-carrier k
            sigma_PkF_GkmF = 0; % Initialize to zero

            % Loop through all femtocells to calculate interference
            for current_f_index = 1:total_femto_count
                % create a 2x2 matrix of current mue coordinates and the
                % coordinates of the interfering femtocell. 
                coord_matrix = [mue_X_coords(mue_index),mue_Y_coords(mue_index);femto_X_coords(current_f_index),femto_Y_coords(current_f_index)];

                % calculate the distance between those two points
                d_femto = pdist(coord_matrix,'euclidean');

                % Path Loss
                PL_femto = 28.0 + 35*log10(d_femto);

                % Channel Gain
                CG_femto = 10^-(PL_femto/10);

                sigma_PkF_GkmF = sigma_PkF_GkmF + (transmitPower_femto*CG_femto);
            end

            % add the femtocell interferers to the denom 
            denominator = denominator + sigma_PkF_GkmF;

            % Loop through macrocells to calculate the gain of each MUE
            for macro_index = 1:Num_Mc
                % calculate the starting and ending indices for the MUE array
                % A(1-150), B(151-300), C(301-450), D(451-600), E(601-750),
                % F(751-900), G(901-1050)
                first_mue_index = (macro_index - 1)*n_mue + 1;
                last_mue_index  = macro_index*n_mue;

                % If this macro_index is desired macro index (mue location is
                % closest to it)
                if (mue_index >= first_mue_index && mue_index <= last_mue_index)
                    % find the distance from every MUE to its host macrocell
                    % create a 2x2 matrix of current MUE coordinates and the
                    % coordinates of the host macrocell
                    coord_matrix = [mue_X_coords(mue_index),mue_Y_coords(mue_index);macro_X_coords(macro_index),macro_Y_coords(macro_index)];

                    % calculate the distance between those two points
                    d_mue = pdist(coord_matrix,'euclidean');

                    % calculate the PL of the MUE based on that distance
                    PL_mue = 28.0 + 35*log10(d_mue);

                    % calculate channel gain for MUE
                    CG_mue = 10^-(PL_mue/10);

                    % numerator of the SINR equation
                    numerator = transmitPower_macro * CG_mue;

                    % combine values into SINR
                    SINR_km = numerator / denominator;

                    % calculate channel capacity
                    channelCapacity_macro = delta_f * log2(1 + (alpha * SINR_km));

                    % multiply beta and the channel capacity
                    % each MUE is assigned to 8 subcarriers, so multiply channel capacity by 8
                    % (because 8 of the Beta values will be 1 and we are summing over all
                    % subcarriers)
                    throughput_macro = channelCapacity_macro * 8;

                    total_throughput = total_throughput + throughput_macro;
                end
            end
        end
    end
        
    % assign the throughput for this femtocell increment to the array
    throughput_macro_array(total_femto_count) = total_throughput;
end

figure; 
plot(femtocell_array,throughput_macro_array, 'o');
xlabel('Number of Femtocells');
ylabel('Throughput (bps)');
title('Teresa plot');

% ^^^ Teresa's Additions ^^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%**************************************************************************
% Femtocell - SINR, Capacity, and Throughput
%**************************************************************************
% From Paper: "After setting the subbands distribution among macrocell regions,
% all the subbands except subband G are reused by femtocells when
% they are deployed inside the macrocell coverage area."


% NOTE: Femtocell distances were already computed in the above
% macrocell section of the code. Those distances are to be reused for this
% section. 

% Precompute Femto User Equipment (FUE) distance vector for simulation
d_FUE_vec  = round(rand(1,f_users)*(macro_radius)+1);

idx=0;
SubCarriers_Assigned = 0;
Tm_femto_vec = zeros(1,length(Num_Fc));
% Loop through all macro cells placing 150 UEs and calculating throughput
% for all 7 macrocells.
for Nm=1:Num_Mc
    
    SubCarriers_Assigned = 0;
    Thruput_FUE_vec = zeros(1,f_users);
    for f=1:f_users
        
        % Convert Macrocell power to watts
        MC_TxP = MC_TxP_vec(1);
        MC_TxP_W = 10^(MC_TxP/10);
        
        % Convert Femtocell power to watts
        FC_TxP_W = 10^(transmitPower_femto/10);
        
%         % Summation of M neighboring Macro-cell's Power & Gain products on sub-carrier k
%         % Equation 4 - Denomonator middle summation
%         sigma_PMp_GMp = 0; % Initialize to zero
%         for mc=1:(Num_Mc-1)
%             d=866;
%             PL_macro = 28.0 + 35*log10(d);
%             Gain = 10^-(PL_macro/10);
%             sigma_PMp_GMp = sigma_PMp_GMp + (MC_TxP_W*Gain);
%         end
        
            
        % Summation of F neighboring Femto-cell Power & Gain products on sub-carrier k
        sigma_PF_GF   = 0; % Initialize to zero
        for fc=1:(Num_Fc)
            PL_femto = 28.0 + 35*log10(d_femto_vec(fc));
            Gain = 10^-(PL_femto/10);
            sigma_PF_GF = sigma_PF_GF + (FC_TxP_W*Gain);
        end
            
         
        % Calculate channel gain (NOTE: Removing the Xsigma and the |H|
        % Rayleigh Gaussian distribution).
        PL_FUE = 28.0 + 35*log10(d_FUE_vec(f));
        Ch_Gain_W = 10^(-PL_FUE/10);
        
        sigma_PMp_GMp = MC_TxP_W * Ch_Gain_W;
        
        % SINR equation for a given Macro-cell on sub-carrier k
        % NOTE: this is equation 4 from the paper
        SINRfk = (FC_TxP_W*Ch_Gain_W)/(10^((Noise_PSD*delta_f)/10) + sigma_PMp_GMp + sigma_PF_GF);
        
        % Capacity of macro user m on sub-carrier k
        Cfk = delta_f*log2(1+alpha*SINRfk);
        

        % Calculate Throughput of the Macro-cell across all m users and all
        % subcarriers available to user.
        Tf = 0; % Initialize to zero
        for k=1:SCpRB
            
            % Create subcarrier assignments, but only assign up to the
            % alotted maximum number based on bandwidth and channel spacing
            if SubCarriers_Assigned < Max_Subcarriers
                Beta_km = 1;
                SubCarriers_Assigned = SubCarriers_Assigned + 1;
            else
                Beta_km = 0;
            end
        
            Tf = Tf + Cfk * Beta_km;
        end
        
        % Assign throughput based on macrouser
        Thruput_FUE_vec(f) = Tf;
        
    end
    
    idx = idx+1;
    Tm_femto_vec(idx) = sum(Thruput_FUE_vec);
    
end

%**************************************************************************
% Figures/Plots
%**************************************************************************
figure; 
plot((1:Num_Mc)*30,(Tm_femto_vec/1e6), 'o');
xlabel('Number of Femto-cells');
ylabel('Throughput (Mbps)');
title('Femtocell Throughput');

