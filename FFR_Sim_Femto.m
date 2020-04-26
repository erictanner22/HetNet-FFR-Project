% Constants
n_femto_per_macro   = 20;      % Number of femtocells 30
n_femto_total       = 140;     % total femtocells (number/macro * 7) 210
r_macro             = 500;     % Radius of Hexagon
n_fue               = 50;      % number of UEs per macrocell 150
number_of_runs      = 100;

throughput_femto_array      = zeros(1,n_femto_total);
throughput_femto_array_avg	= zeros(1,n_femto_total);
femtocell_array             = 1:n_femto_total;

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

%choose only N_femto points for femtocells
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

% combine all femto coordinates from every macrocell into one array
femto_X_coords = [Ac_x Bc_x Cc_x Dc_x Ec_x Fc_x Gc_x];
femto_Y_coords = [Ac_y Bc_y Cc_y Dc_y Ec_y Fc_y Gc_y];

% increment total femtocells for graphing
% Only include femtocells within macrocell A
for total_femto_count = 1:n_femto_total

    % hold the sum of all MUE throughput - reset at every femtocell
    % increment
    total_throughput = 0;
   
    % run Monte Carlo Simulation 100 times at every increment of femtocell
    for i = 1:number_of_runs
    
        %--------- Generate new locations for MUEs ----------------------%
        %The method used here is to generate many points in a square and choose N points that fall within the hexagon
        %Generate 3*n_mue random points with square that is 2R by 2R
        A_fue_x = (r_macro-rand(1, 3*n_fue)*2*r_macro) + A_center_X;
        A_fue_y = (r_macro-rand(1, 3*n_fue)*2*r_macro) + A_center_Y;
        B_fue_x = (r_macro-rand(1, 3*n_fue)*2*r_macro) + B_center_X;
        B_fue_y = (r_macro-rand(1, 3*n_fue)*2*r_macro) + B_center_Y;
        C_fue_x = (r_macro-rand(1, 3*n_fue)*2*r_macro) + C_center_X;
        C_fue_y = (r_macro-rand(1, 3*n_fue)*2*r_macro) + C_center_Y;
        D_fue_x = (r_macro-rand(1, 3*n_fue)*2*r_macro) + D_center_X;
        D_fue_y = (r_macro-rand(1, 3*n_fue)*2*r_macro) + D_center_Y;
        E_fue_x = (r_macro-rand(1, 3*n_fue)*2*r_macro) + E_center_X;
        E_fue_y = (r_macro-rand(1, 3*n_fue)*2*r_macro) + E_center_Y;
        F_fue_x = (r_macro-rand(1, 3*n_fue)*2*r_macro) + F_center_X;
        F_fue_y = (r_macro-rand(1, 3*n_fue)*2*r_macro) + F_center_Y;
        G_fue_x = (r_macro-rand(1, 3*n_fue)*2*r_macro) + G_center_X;
        G_fue_y = (r_macro-rand(1, 3*n_fue)*2*r_macro) + G_center_Y;

        %There is a command in MATLAB inploygon. 
        %The command finds points within a polygon region.
        %get the MUE points within the polygon
        A_IN_fue = inpolygon(A_fue_x, A_fue_y, Av_x, Av_y);
        B_IN_fue = inpolygon(B_fue_x, B_fue_y, Bv_x, Bv_y);
        C_IN_fue = inpolygon(C_fue_x, C_fue_y, Cv_x, Cv_y);
        D_IN_fue = inpolygon(D_fue_x, D_fue_y, Dv_x, Dv_y);
        E_IN_fue = inpolygon(E_fue_x, E_fue_y, Ev_x, Ev_y);
        F_IN_fue = inpolygon(F_fue_x, F_fue_y, Fv_x, Fv_y);
        G_IN_fue = inpolygon(G_fue_x, G_fue_y, Gv_x, Gv_y);

        %drop MUE nodes outside the hexagon
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

        %choose only n_mue points for MUEs
        A_idx = randperm(length(A_fue_x));
        A_fue_x = A_fue_x(A_idx(1:n_fue));
        A_fue_y = A_fue_y(A_idx(1:n_fue));
        B_idx = randperm(length(B_fue_x));
        B_fue_x = B_fue_x(B_idx(1:n_fue));
        B_fue_y = B_fue_y(B_idx(1:n_fue));
        C_idx = randperm(length(C_fue_x));
        C_fue_x = C_fue_x(C_idx(1:n_fue));
        C_fue_y = C_fue_y(C_idx(1:n_fue));
        D_idx = randperm(length(D_fue_x));
        D_fue_x = D_fue_x(D_idx(1:n_fue));
        D_fue_y = D_fue_y(D_idx(1:n_fue));
        E_idx = randperm(length(E_fue_x));
        E_fue_x = E_fue_x(E_idx(1:n_fue));
        E_fue_y = E_fue_y(E_idx(1:n_fue));
        F_idx = randperm(length(F_fue_x));
        F_fue_x = F_fue_x(F_idx(1:n_fue));
        F_fue_y = F_fue_y(F_idx(1:n_fue));
        G_idx = randperm(length(G_fue_x));
        G_fue_x = G_fue_x(G_idx(1:n_fue));
        G_fue_y = G_fue_y(G_idx(1:n_fue));
        
        % combine all MUE coordinates from every macrocell into one array
        fue_X_coords = [A_fue_x B_fue_x C_fue_x D_fue_x E_fue_x F_fue_x G_fue_x];
        fue_Y_coords = [A_fue_y B_fue_y C_fue_y D_fue_y E_fue_y F_fue_y G_fue_y];
        %----------------------------------------------------------------%        
        
        % Loop through all FUEs
        for fue_index = 1:n_fue
            
            %-----------Find the nearest femtocell to assign the FUE to------%
            % initialize shortest distance as the diameter of the network
            d_shortest = 3000;
            assigned_femtocell_index = 500;
            
            % loop through all femtocells in this iteration (1-210)
            for femto_index = 1:total_femto_count
                % create a 2x2 matrix of current mue coordinates and the
                % coordinates of the interfering macrocell. 
                coord_matrix = [fue_X_coords(fue_index),fue_Y_coords(fue_index);femto_X_coords(femto_index),femto_Y_coords(femto_index)];

                % calculate the distance between those two points
                d_femto = pdist(coord_matrix,'euclidean');

                % if the current distance is less than the saved max distance
                if d_femto <= d_shortest
                    % replace the new shortest distance
                    d_shortest = d_femto;

                    assigned_femtocell_index = femto_index;
                end
            end
            %----------------------------------------------------------------%
            
            % multiply noise spectral density by the subcarrier spacing and convert
            % to Watts
            denominator = 10^((Noise_PSD * delta_f)/10);

            % calculate macrocell interference
            % Summation of M neighboring Macro-cell's Power & Gain products on sub-carrier k
            sigma_Pkm_GkfM = 0; % Initialize to zero

            % Calculate interference to MUE from interferer macrocells
            for current_macro_index = 1:Num_Mc
                % create a 2x2 matrix of current mue coordinates and the
                % coordinates of the interfering macrocell. 
                coord_matrix = [fue_X_coords(fue_index),fue_Y_coords(fue_index);macro_X_coords(current_macro_index),macro_Y_coords(current_macro_index)];

                % calculate the distance between those two points
                d_macro = pdist(coord_matrix,'euclidean');

                % outdoor pathloss - equation (2) from paper
                PL_macro = 28.0 + 35*log10(d_macro);

                % equation (3) from paper
                CG_macro = 10^-(PL_macro/10);

                % Add up all the interferers
                sigma_Pkm_GkfM = sigma_Pkm_GkfM + (transmitPower_macro*CG_macro);
            end

            % add the macrocell interferers to the denom 
            denominator = denominator + sigma_Pkm_GkfM;

            % femtocell interference
            % Summation of F neighboring Femto-cell Power & Gain products on sub-carrier k
            sigma_PkF_GkfF = 0; % Initialize to zero

            % Loop through all femtocells to calculate interference
            for current_f_index = 1:total_femto_count
                
                % don't calculate the interference of the assigned f_cell
                if not(current_f_index == assigned_femtocell_index)
                    % create a 2x2 matrix of current mue coordinates and the
                    % coordinates of the interfering femtocell. 
                    coord_matrix = [fue_X_coords(fue_index),fue_Y_coords(fue_index);femto_X_coords(current_f_index),femto_Y_coords(current_f_index)];

                    % calculate the distance between those two points
                    d_femto = pdist(coord_matrix,'euclidean');

                    % Path Loss
                    PL_femto = 28.0 + 35*log10(d_femto);

                    % Channel Gain
                    CG_femto = 10^-(PL_femto/10);

                    sigma_PkF_GkfF = sigma_PkF_GkfF + (transmitPower_femto*CG_femto);
                end
            end

            % add the femtocell interferers to the denom 
            denominator = denominator + sigma_PkF_GkfF;

            %--------Numerator--------------------------------------%
            
            % find the distance from every MUE to its host macrocell
            % create a 2x2 matrix of current MUE coordinates and the
            % coordinates of the host macrocell
            coord_matrix = [fue_X_coords(fue_index),fue_Y_coords(fue_index);femto_X_coords(assigned_femtocell_index),femto_Y_coords(assigned_femtocell_index)];

            % calculate the distance between those two points
            d_fue = pdist(coord_matrix,'euclidean');

            % calculate the PL of the MUE based on that distance
            PL_fue = 28.0 + 35*log10(d_fue);

            % calculate channel gain for MUE
            CG_fue = 10^-(PL_fue/10);

            % numerator of the SINR equation
            numerator = transmitPower_femto * CG_fue;

            % combine values into SINR
            SINR_kf = numerator / denominator;

            % calculate channel capacity
            channelCapacity_femto = delta_f * log2(1 + (alpha * SINR_kf));

            % multiply beta and the channel capacity
            % each MUE is assigned to 8 subcarriers, so multiply channel capacity by 8
            % (because 8 of the Beta values will be 1 and we are summing over all
            % subcarriers)
            throughput_femto = channelCapacity_femto * 8;

            total_throughput = total_throughput + throughput_femto;
        end
    end
        
    % assign the throughput for this femtocell increment to the array
    throughput_femto_array(total_femto_count) = total_throughput;
    
    throughput_femto_array_avg(total_femto_count) = (total_throughput/number_of_runs);
end


figure; 
plot(femtocell_array,throughput_femto_array_avg, 'x');
xlabel('Number of Femtocells');
ylabel('Averaged Throughput (bps)');
title('Throughput of the UEs connected with Femtocell');


figure; 
plot(femtocell_array,throughput_femto_array, 'x');
xlabel('Number of Femtocells');
ylabel('Throughput (bps)');
title('Throughput of the UEs connected with Femtocell');

% create an array to plot every 10 femto increments. 
averaged_throughput_array = zeros(1,(n_femto_total/10));