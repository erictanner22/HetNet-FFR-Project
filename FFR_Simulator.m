%**************************************************************************
% Filename: FFR_Simulator.m
% Group Name: GTW-E
% Date: 03/17/2020
% Description: Script used to simulate Fractional Frequency Reuse(FFR)
% medthods for HetNet applications. Four FFR algorithms are modeled to
% analyze performance tradeoffs between algorithm approaches. 
%
% FFR Algorithms:
%                FFR-3SL (Proposed Paper) 
%                FFR-3  (Reference:10) 
%                OSFFR  (Reference:11)
%                FFR-3R (Reference:12) 
%
%**************************************************************************

%----------------------------------------------------
% Common Sim Variables
%----------------------------------------------------
BER     = 10^-6;         % Target Bit Error Rate(BER)
alpha   = -1.5/log(BER); % Constant for target BER
deltaf  = 15e3;          % Subcarrier spacing (Hz)
CH_BW   = 20e6;          % Channel Bandwidth (Hz)
m_users = 150;           % Macrocell Users
f_users = 150;           % Femtocell Users
SCpRB   = 12;            % Sub-carriers per resource block
Num_RB  = 100;           % Total number of resource blocks
Num_SC  = Num_RB*SCpRB;  % Number of subcarriers
MC_TxP  = [10 15 20];    % Macrocell Base Station Transmit Power
FC_TxP  = 20e-3;         % Femtocell Base Station Transmit Power
No_PSD  = -174;          % Noise Power Spectral Density (dBm/Hz)
Nm      = 7;             % Number of Macrocells
Nf      = 210;           % Number of Femtocells
d_vec  = 5;              % Distance from base station to user in meters
Lwalls = [7 10 15];      % Loss through walls [light internal, internal, external]
wall_type = 1;           % Selects wall type from array of wall loss vector
                         % hard coded to 1(light internal) for now, will 
                         % implement selector code later.

                         

PL_outdoor = 28.0 + 35*log10(d_vec);
PL_indoor  = 38.5 + 20*log10(d_vec)+Lwalls(wall_type);
PL_vec = [PL_outdoor PL_indoor];
PL_type = 1;        % Selects PL type being either indoor or outdoor.
                    % hard coded to 1(outdoor) for now, will implement selector code
                    % later.
PL = PL_vec(PL_type);
                    
%--------------------------------------------------------------------------
% FFR-3SL Code (Proposed Paper)
%--------------------------------------------------------------------------






%--------------------------------------------------------------------------
% FFR-3 Code FFR-3  (Reference:10)
%--------------------------------------------------------------------------


%-------------------------------------------
% Macrocell - SINR, Capacity, and Throughput
%-------------------------------------------
idx = 0;
Nm_vec = 5:5:200;
for Num_Mc = Nm_vec
    
    Ch_Gain = 10^-(PL/10);
    
    MC_TxP = MC_TxP(1);
    MC_TxP_W = 10^(MC_TxP/10);
    
    FC_TxP_W = 10^(FC_TxP/10);
    
    % Summation of M neighboring Macro-cell's Power & Gain products on sub-carrier k
    sigma_PMp_GMp = 0; % Initialize to zero
    for m=1:(Num_Mc-1)
        sigma_PMp_GMp = sigma_PMp_GMp + (MC_TxP_W*Ch_Gain);
    end
    
    % Summation of F neighboring Femto-cell Power & Gain products on sub-carrier k
    sigma_PF_GF = 0; % Initialize to zero
    %for f=1:(Nf-1)
    for f=1:(Num_Mc)
        sigma_PF_GF = sigma_PF_GF + (FC_TxP_W*Ch_Gain);
    end
    
    % SINR equation for a given Macro-cell on sub-carrier k
    SINRmk = (MC_TxP_W*Ch_Gain)/(10^((No_PSD*deltaf)/10) + sigma_PMp_GMp + sigma_PF_GF);
    
    % Capacity of macro user m on sub-carrier k
    Cmk = deltaf*log2(1+alpha*SINRmk);
    
    % Calculate Throughput of the Macro-cell across all m users
    Tm = 0; % Initialize to zero
    for m=1:m_users
        Tm = Tm + Cmk;
    end
    
    % Store the throughput into the vector (vector used for plotting)
    idx = idx+1;
    Tm_vec(idx) = Tm;
    
end

%-------------------------------------------
% Femtocell - SINR, Capacity, and Throughput
%-------------------------------------------
% idx = 0;
% Nm_vec = 5:5:200;
% for Num_Mc = Nm_vec
%     
%     idx = idx+1;
%     
%     Ch_Gain = 10^-(PL/10);
%     
%     MC_TxP = MC_TxP(1);
%     MC_TxP_W = 10^(MC_TxP/10);
%     
%     FC_TxP_W = 10^(FC_TxP/10);
%     
%     % Summation of M neighboring Macro-cell's Power & Gain products on sub-carrier k
%     sigma_PMp_GMp = 0; % Initialize to zero
%     for m=1:(Num_Mc-1)
%         sigma_PMp_GMp = sigma_PMp_GMp + (MC_TxP_W*Ch_Gain);
%     end
%     
%     % Summation of F neighboring Femto-cell Power & Gain products on sub-carrier k
%     sigma_PF_GF = 0; % Initialize to zero
%     for f=1:(Nf-1)
%         sigma_PF_GF = sigma_PF_GF + (FC_TxP_W*Ch_Gain);
%     end
%     
%     % NOTE: Needs more work.
%     % % SINR equation for a given Femto-cell on sub-carrier k
%     % SINRfk = (FC_TxP_W*Ch_Gain)/(10^((No_PSD*deltaf)/10) + sigma_PMp_GMp + sigma_PF_GF);
%     
% 
% end


%--------------------------------------------------------------------------
% OSFFR Code (Reference:11)
%--------------------------------------------------------------------------








%--------------------------------------------------------------------------
% FFR-3R Code FFR-3R (Reference:12)
%--------------------------------------------------------------------------






%--------------------------------------------------------------------------
% Analysis
%--------------------------------------------------------------------------





%--------------------------------------------------------------------------
% Figures/Plots
%--------------------------------------------------------------------------

figure; 
plot(Nm_vec,Tm_vec, 'o');
xlabel('Number of Macro-cells');
ylabel('Throughput (bps)');
