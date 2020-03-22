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
BER    = 10^-6;         % Target Bit Error Rate(BER)
alpha  = -1.5/log(BER); % Constant for target BER
deltaf = 15e3;          % Subcarrier spacing (Hz)
CH_BW  = 5e6;           % Channel Bandwidth (Hz)
M      = 0;             % Macro cell M 
m      = 0;             % Macro user m
k      = 0;             % Subcarrier k
F      = 0;             % Femto cell F
f      = 0;             % Femto user f

%--------------------------------------------------------------------------
% FFR-3SL Code (Proposed Paper)
%--------------------------------------------------------------------------






%--------------------------------------------------------------------------
% FFR-3 Code FFR-3  (Reference:10) 
%--------------------------------------------------------------------------
d_vec  = (1:5:100); % Distance from base station to user in meters
Lwalls = [7 10 15]; % Loss through walls [light internal, internal, external]
wall_type = 1;      % Selects wall type from array of wall loss vector
                    % hard coded to 1(light internal) for now, will 
                    % implement selector code later.
                    
PL_outdoor = 28.0 + 35*log10(d);
PL_indoor  = 38.5 + 20*log10(d)+Lwalls(wall_type);
PL_vec = [PL_outdoor PL_indoor];

PL_type = 1;        % Selects PL type being either indoor or outdoor.
                    % hard coded to 1(outdoor) for now, will implement selector code
                    % later.

Ch_Gain = 10^(-PL_vec(PL_type)/10);




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

