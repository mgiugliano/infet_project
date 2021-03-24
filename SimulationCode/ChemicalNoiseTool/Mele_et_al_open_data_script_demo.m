%   Demo source file from article entitled:
%   "General Model and Equivalent Circuit for the Chemical Noise Spectrum 
%    Associated to Surface Charge Fluctuation in Potentiometric Sensors"
%
%   Journal: IEEE Sensors Journal
%
%   DOI: 10.1109/JSEN.2020.3038036
%   
%   Authors: L. J. Mele, P. Palestri, L. Selmi
%
%   Date: 19/02/2021

clear all; close all; clc

% Constants
q      = 1.602176634e-19;      % C
k_B    = 1.380649e-23;         % J/K
Temp   = 298.16;               % K

% Definition of symbol s for the s-polynomial Y_S
s=sym('s');

% Plot settings
linestyle = {'-';':';'--'};

% Frequency vector
freq = logspace(-6,6,500);

% Initializating the surface charge density PSD vector
S_QsQs = zeros(1,length(freq));

% pH values
pH = [5, 7, 9];

% Parameter values taken from Zhang et al. 2018, ref. [21]
Ns = 9e18;
Ka = 2e-8;
Kb = 2e-5;
ka_f = 6.9554e+07;  % @f=10Hz
kb_f = 267.5383;
WL = 7.854e-5;     % Sensing layer area [m^2]

for index_i = 1:length(pH)
    
    % pH selection, data from Zhang et al. 2018
    switch pH(index_i)
        case 5
            Hs = 0.81e-6;
        case 7
            Hs = 9e-6;
        case 9
            Hs = 0.4e-6;
    end
    
    %% Computation of the Surface charge density PSD, (Eq. 14)
    
    % f0_ denotes f^{'0}, the truncated vector of state ocupation functions
    % at equilibrium (f_{MOH} is the eliminated state)
    f0_ =  [(Ka*Kb)/(Hs^2 + Kb*Hs + Ka*Kb); ...
               Hs^2/(Hs^2 + Kb*Hs + Ka*Kb)];
    
    % Setting the matrix G, T, R, Omega
    G = [ -Hs*ka_f,             Ka*ka_f,        0;
           Hs*ka_f, - Hs*kb_f - Ka*ka_f,  Kb*kb_f;
                 0,             Hs*kb_f, -Kb*kb_f];
    T = [1     0     0;  0     0     1];
    R = [1     0 ; -1    -1; 0     1];
    Omega   = T*G*R;
    
    % Truncated vector of states charge
    z_  = [-1 1];
    
    % Gamma matrix calculus (Eq. 15)
    Gamma = f0_.*eye(2)-f0_*f0_';
    
    % Function implementing (sI - Omega)
    V = @(s) eye(2).*s-Omega;
    
    for freq_index = 1:length(freq)
        jw = 2*pi*freq(freq_index)*1i;
        S_QsQs(freq_index) = 4*q^2*Ns/WL*z_*(Gamma*real(inv(V(jw)))')*z_';
    end
    
    % Plot
    figure(1)
    p(index_i) = loglog(freq,S_QsQs,'-','LineStyle',linestyle{index_i});
    hold on
    
    
    %% Compuation of the equivalent admittance Y_S, according to Eq. 29
    % Example shown for the case pH = 7
    if pH(index_i) == 7
        I = eye(2);
        
        G_1 = [ -ka_f,     0,    0;
                 ka_f,     0,    0;
                    0,     0,    0];
        
        G_2 = [     0,     0,    0;
                    0, -kb_f,    0;
                    0,  kb_f,    0];
        
        Zeta_1 = [ 1,  0,  0;
                   1,  0,  0;
                   0,  0,  0];
        
        Zeta_2 = [ 0,  0,  0;
                   0,  1,  0;
                   0,  1,  0];
        
        A_1 = [   Hs,  0,  0;
                  Hs,  0,  0;
                   0,  0,  0];
        
        A_2 = [    0,  0,  0;
                   0, Hs,  0;
                   0, Hs,  0];
        
        % Full vector of states' occupation probability at equilibrium
        f0 =  [(Ka*Kb)/(Hs^2 + Kb*Hs + Ka*Kb);
               (Hs*Kb)/(Hs^2 + Kb*Hs + Ka*Kb);
                  Hs^2/(Hs^2 + Kb*Hs + Ka*Kb)];
        
        Theta = [1, -1; -1, 1];
        
        % Equation 29
        Y_S_vector = q^2*Ns/k_B/Temp*((I-(Theta.*Omega)/s)\(z_'.*( T*( ...
            G_1*Zeta_1*A_1 + G_2*Zeta_2*A_2)*f0)));
        
        disp('Eq. 30 in the manuscript, calculated at pH = 7:')
        % Print Y_S as a Laplace polynomial
        Y_S = ones(1,2)*Y_S_vector
    end
end

% Editing Axes
grid on
xlim([1e-6,1e6])
ylim([1e-27,1e-12])
set(gca, 'XTick', 10.^(-6:2:6))
set(gca, 'YTick', 10.^(-27:3:-12))

% Title and Labels
set(gca,'XMinorTick','off','YMinorTick','off')
xlabel('f [Hz]');
ylabel('S_{Q_sQ_s} [C^2/m^4/Hz]');

% Set legend
legend(p,strcat('pH = ',num2str([9 7 5]')));