% Constants
noise_PSD           = -174;             % Noise Power Spectral Density (dBm/Hz)
BER                 = 10^-4;            % Target Bit Error Rate(BER)
alpha               = -1.5/log(BER);	% Constant for target BER
delta_f             = 15e3;             % Subcarrier spacing (Hz)
transmitPower_macro	= 20;               % Macrocell transmit power in Watts
transmitPower_femto = 20e-3;            % Femtocell Base Station Transmit Power in Watts
CH_BW               = 20e6;             % Channel Bandwidth (Hz)
max_subcarriers     = CH_BW/delta_f;    % Maximum number of subcarriers

n_macro_total       = 7;        % Number of Macrocells
n_femto_per_macro   = 30;       % Number of femtocells 30
n_femto_total       = 210;      % total femtocells (number/macro * 7) 210
r_macro             = 500;      % Radius of Hexagon
n_fue_per_macro     = 150;      % number of UEs per macrocell 150
n_fue_total         = 1050;     % total femtocells in the network
number_of_runs      = 100;      % number of times to run the Monte Carlo sim

% Store all possible increments of the femtocell count to increment through
femtocell_array     = 1:n_femto_total;

% Initialize the throughput arrays
throughput_femto_array      = zeros(1,n_femto_total);
throughput_femto_array_avg	= zeros(1,n_femto_total);

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
for total_femto_count = 1:n_femto_total

    % Hold the sum of all FUE throughput - reset at every femtocell
    % increment
    total_throughput = 0;
    
    % Initialize the subcarrier assignment counter to 0
    subcarriers_assigned = 0;
    
    % Run Monte Carlo Simulation 100 times at every increment of femtocell
    for i = 1:number_of_runs
    
        %--------- Generate new locations for FUEs ----------------------%
        % Generate many points in a square and choose N points that fall within the hexagon
        % Generate 3*n_fue_per_macro random points with square that is 2R by 2R
        A_fue_x = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + A_center_X;
        A_fue_y = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + A_center_Y;
        B_fue_x = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + B_center_X;
        B_fue_y = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + B_center_Y;
        C_fue_x = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + C_center_X;
        C_fue_y = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + C_center_Y;
        D_fue_x = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + D_center_X;
        D_fue_y = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + D_center_Y;
        E_fue_x = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + E_center_X;
        E_fue_y = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + E_center_Y;
        F_fue_x = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + F_center_X;
        F_fue_y = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + F_center_Y;
        G_fue_x = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + G_center_X;
        G_fue_y = (r_macro-rand(1, 3*n_fue_per_macro)*2*r_macro) + G_center_Y;

        % The MATLAB command inploygon finds points within a polygon region
        A_IN_fue = inpolygon(A_fue_x, A_fue_y, Av_x, Av_y);
        B_IN_fue = inpolygon(B_fue_x, B_fue_y, Bv_x, Bv_y);
        C_IN_fue = inpolygon(C_fue_x, C_fue_y, Cv_x, Cv_y);
        D_IN_fue = inpolygon(D_fue_x, D_fue_y, Dv_x, Dv_y);
        E_IN_fue = inpolygon(E_fue_x, E_fue_y, Ev_x, Ev_y);
        F_IN_fue = inpolygon(F_fue_x, F_fue_y, Fv_x, Fv_y);
        G_IN_fue = inpolygon(G_fue_x, G_fue_y, Gv_x, Gv_y);

        % Drop FUE nodes outside the hexagon
        A_fue_x = A_fue_x(A_IN_fue);
        A_fue_y = A_fue_y(A_IN_fue);
        B_fue_x = B_fue_x(B_IN_fue);
        B_fue_y = B_fue_y(B_IN_fue);
        C_fue_x = C_fue_x(C_IN_fue);
        C_fue_y = C_fue_y(C_IN_fue);
        D_fue_x = D_fue_x(D_IN_fue);
        D_fue_y = D_fue_y(D_IN_fue);
        E_fue_x = E_fue_x(E_IN_fue);
        E_fue_y = E_fue_y(E_IN_fue);
        F_fue_x = F_fue_x(F_IN_fue);
        F_fue_y = F_fue_y(F_IN_fue);
        G_fue_x = G_fue_x(G_IN_fue);
        G_fue_y = G_fue_y(G_IN_fue);

        % Choose only n_fue_per_macro points for FUEs
        A_idx = randperm(length(A_fue_x));
        A_fue_x = A_fue_x(A_idx(1:n_fue_per_macro));
        A_fue_y = A_fue_y(A_idx(1:n_fue_per_macro));
        B_idx = randperm(length(B_fue_x));
        B_fue_x = B_fue_x(B_idx(1:n_fue_per_macro));
        B_fue_y = B_fue_y(B_idx(1:n_fue_per_macro));
        C_idx = randperm(length(C_fue_x));
        C_fue_x = C_fue_x(C_idx(1:n_fue_per_macro));
        C_fue_y = C_fue_y(C_idx(1:n_fue_per_macro));
        D_idx = randperm(length(D_fue_x));
        D_fue_x = D_fue_x(D_idx(1:n_fue_per_macro));
        D_fue_y = D_fue_y(D_idx(1:n_fue_per_macro));
        E_idx = randperm(length(E_fue_x));
        E_fue_x = E_fue_x(E_idx(1:n_fue_per_macro));
        E_fue_y = E_fue_y(E_idx(1:n_fue_per_macro));
        F_idx = randperm(length(F_fue_x));
        F_fue_x = F_fue_x(F_idx(1:n_fue_per_macro));
        F_fue_y = F_fue_y(F_idx(1:n_fue_per_macro));
        G_idx = randperm(length(G_fue_x));
        G_fue_x = G_fue_x(G_idx(1:n_fue_per_macro));
        G_fue_y = G_fue_y(G_idx(1:n_fue_per_macro));
        
        % Combine all FUE coordinates from every macrocell into one array
        fue_X_coords = [A_fue_x B_fue_x C_fue_x D_fue_x E_fue_x F_fue_x G_fue_x];
        fue_Y_coords = [A_fue_y B_fue_y C_fue_y D_fue_y E_fue_y F_fue_y G_fue_y];
        %----------------------------------------------------------------%        
        
        % Loop through all FUEs in the network
        for fue_index = 1:n_fue_total
            
            %-----------Find the nearest femtocell to assign the FUE to------%
            % Initialize shortest distance as the diameter of the network
            d_shortest_femto = 3000;
            assigned_femtocell_index = 500;
            
            % Loop through all femtocells in this iteration
            for femto_index = 1:total_femto_count
                % create a 2x2 matrix of current FUE coordinates and the
                % coordinates of the current femtocell index
                coord_matrix = [fue_X_coords(fue_index),fue_Y_coords(fue_index);femto_X_coords(femto_index),femto_Y_coords(femto_index)];

                % calculate the distance between those two points
                d_femto = pdist(coord_matrix,'euclidean');

                % if the current distance is less than the saved min distance
                if d_femto <= d_shortest_femto
                    % replace the new shortest distance
                    d_shortest_femto = d_femto;

                    % save the index of the closest femtocell for later use
                    assigned_femtocell_index = femto_index;               
                end
            end
            %----------------------------------------------------------------%
            
            % Start the calculation of Equation (5)
            % multiply noise spectral density by the subcarrier spacing and convert
            % to Watts
            denominator = 10^((noise_PSD * delta_f)/10);

            % calculate macrocell interference
            % Summation of neighboring macrocell's power & gain products on subcarrier k
            sigma_Pkm_GkfM = 0; % Initialize to zero
            
            % Calculate interference to FUE from interferer macrocells
            for current_macro_index = 1:n_macro_total
                % create a 2x2 matrix of current FUE coordinates and the
                % coordinates of the interfering macrocell. 
                coord_matrix = [fue_X_coords(fue_index),fue_Y_coords(fue_index);macro_X_coords(current_macro_index),macro_Y_coords(current_macro_index)];

                % calculate the distance between those two points
                d_macro = pdist(coord_matrix,'euclidean');
                
                % Outdoor Pathloss - Equation (2) from paper
                PL_macro = 28.0 + 35*log10(d_macro);

                % Channel Gain - Equation (3) from paper
                CG_macro = 10^-(PL_macro/10);

                % Add up all the interferers
                sigma_Pkm_GkfM = sigma_Pkm_GkfM + (transmitPower_macro*CG_macro);
            end
            
            % add the macrocell interferers to the denominator of Equation
            % (5)
            denominator = denominator + sigma_Pkm_GkfM;

            % Femtocell interference
            % Summation of F neighboring femtocell power & gain products on sub-carrier k
            sigma_PkF_GkfF = 0; % Initialize to zero

            % Loop through all femtocells to calculate interference
            for current_f_index = 1:total_femto_count
                
                % Don't calculate the interference of the assigned f_cell
                if not(current_f_index == assigned_femtocell_index)
                    % Create a 2x2 matrix of current FUE coordinates and the
                    % Coordinates of the interfering femtocell. 
                    coord_matrix = [fue_X_coords(fue_index),fue_Y_coords(fue_index);femto_X_coords(current_f_index),femto_Y_coords(current_f_index)];

                    % Calculate the distance between those two points
                    d_femto = pdist(coord_matrix,'euclidean');

                    % Outdoor Pathloss - Equation (2) from paper
                    PL_femto = 28.0 + 35*log10(d_femto);

                    % Channel Gain - Equation (3) from paper
                    CG_femto = 10^-(PL_femto/10);

                    sigma_PkF_GkfF = sigma_PkF_GkfF + (transmitPower_femto*CG_femto);
                end
            end

            % Add the femtocell interferers to the denominator of Equation
            % (5)
            denominator = denominator + sigma_PkF_GkfF;

            %--------Numerator of Equation (5)---------------------------%
            
            % Create a 2x2 matrix of current FUE coordinates and the
            % coordinates of the host femtocell
            coord_matrix = [fue_X_coords(fue_index),fue_Y_coords(fue_index);femto_X_coords(assigned_femtocell_index),femto_Y_coords(assigned_femtocell_index)];

            % Calculate the distance between those two points
            d_fue = pdist(coord_matrix,'euclidean');

            % Outdoor Pathloss - Equation (2) from paper
            PL_fue = 28.0 + 35*log10(d_fue);

            % Channel Gain - Equation (3) from paper
            CG_fue = 10^-(PL_fue/10);

            % Numerator of the SINR equation - Equation (5)
            numerator = transmitPower_femto * CG_fue;

            % Combine values into SINR to complete Equation (5)
            SINR_kf = numerator / denominator;

            % Calculate channel capacity - Equation (7) from paper
            channelCapacity_femto = delta_f * log2(1 + (alpha * SINR_kf));

            % Create subcarrier assignments, but only assign up to the
            % alotted maximum number based on bandwidth and channel spacing
            % Equation (9) from paper
            if subcarriers_assigned < (max_subcarriers - 12)
                throughput_femto = channelCapacity_femto * 12;
                subcarriers_assigned = subcarriers_assigned + 12;
            else
                throughput_femto = 0;
            end

            % Summation of the throughput over all FUEs and
            % subcarriers
            total_throughput = total_throughput + throughput_femto;
        end
    end
        
    % Assign the throughput for this femtocell increment to the array
    throughput_femto_array(total_femto_count) = total_throughput;
    
    throughput_femto_array_avg(total_femto_count) = (total_throughput/number_of_runs);
end


figure; 
plot(femtocell_array,throughput_femto_array_avg/1e6, 'x');
xlabel('Number of Femtocells');
ylabel('Averaged Throughput over Runs (Mbps)');
title('Throughput of the UEs connected with Femtocell');


figure; 
plot(femtocell_array,throughput_femto_array/1e6, 'x');
xlabel('Number of Femtocells');
ylabel('Throughput (Mbps)');
title('Throughput of the UEs connected with Femtocell');

figure; 
plot(femtocell_array(1: 10 : end),throughput_femto_array(1: 10 : end)/1e6, 'x');
xlabel('Number of Femtocells');
ylabel('Throughput (Mbps)');
title('Throughput of the UEs connected with Femtocell');