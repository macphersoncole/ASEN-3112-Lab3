%% ASEN 3112 LAB 3 4-Element Model

% Housekeeping
clear;
close all;
clc;

%% Constants

% Longitudinal Dimensions
fuse_span = 22; % Fuselage Span [in]
wing_span = 18; % Wing Span [in]
L = 12; % End of Shaker to Start of Tail Span [in]
L_E = 4.5; % Elevator Span [in]
L_R = 5; % Rudder Span [in]

% Beam Dimansions
w = 1; % Width (every member) [in]
h = 1/8; % Thickness: fuselage & wing [in]
h_E = 1/4; % Thickness: elevator [in]
h_R = 1/25; % Thickness: rudder [in]

% Material Properties (Al 6063-T7 stock)
E = 10.175*(10^6); % Elastic modulus [psi]
rho = 2.505*(10^(-4)); % Density [lb-sec^2/in^4]

%% Calculated Constants

A = w*h; % Cross-sectional area of the fuselage [in^2]
A_E = w*h_E; % Cross-sectional area of the elevator [in^2]
A_G = w*h_R; % Cross-sectional area of the rudder [in^2]

c_M4 = (rho*A*L)/(806400); % M matrix coefficient
I_ZZ = (w*(h^3))/(12); % Moment of Inertia
c_K4 = (8*E*I_ZZ)/(L^3); % K matrix coefficient

M_T = 1.131*rho; % Mass of tail assembly
S_T = 0.5655*rho; % First mass-moment of tail assembly wrt B'
I_T = 23.124*rho; % Second mass-moment of tail assembly wrt B'

M_4 = (c_M4)*[77088 2916*L 23712 -1284*L 0 0 0 0 0 0;...
    2916*L 172*(L^2) 1284*L -73*(L^2) 0 0 0 0 0 0;...
    23712 1284*L 154176 0 23712 -1284*L 0 0 0 0;...
    -1284*L -73*(L^2) 0 344*(L^2) 1284*L -73*(L^2) 0 0 0 0;...
    0 0 23712 1284*L 154176 0 23712 -1284*L 0 0;...
    0 0 -1284*L -73*(L^2) 0 344*(L^2) 1284*L -73*(L^2) 0 0;...
    0 0 0 0 23712 1284*L 154176 0 23712 -1284*L;...
    0 0 0 0 -1284*L -73*(L^2) 0 344*(L^2) 1284*L -73*(L^2);...
    0 0 0 0 0 0 27312 1284*L 77088 -2916*L;...
    0 0 0 0 0 0 -1284*L -73*(L^2) -2916*L 172*(L^2)] + ...
    [0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 M_T S_T;...
    0 0 0 0 0 0 0 0 S_T I_T]; % M matrix

K_4 = (c_K4)*[96 12*L -96 12*L 0 0 0 0 0 0;...
    12*L 2*(L^2) -12*L L^2 0 0 0 0 0 0;...
    -96 -12*L 192 0 -96 12*L 0 0 0 0;...
    12*L L^2 0 4*(L^2) -12*L L^2 0 0 0 0;...
    0 0 -96 -12*L 192 0 -96 12*L 0 0;...
    0 0 12*L L^2 0 4*(L^2) -12*L L^2 0 0;...
    0 0 0 0 -96 -12*L 192 0 -96 12*L;...
    0 0 0 0 12*L L^2 0 4*(L^2) -12*L L^2;...
    0 0 0 0 0 0 -96 -12*L 96 -12*L;...
    0 0 0 0 0 0 12*L L^2 -12*L 2*(L^2)]; % K matrix

M_hat4 = M_4(3:10,3:10); % Reduced M matrix
K_hat4 = K_4(3:10,3:10); % Reduced K matrix

[U,omega_sq] = eig(K_hat4,M_hat4); % Eigen calculation

omega = sqrt(omega_sq); % Finding omega

freq = (omega)/(2*pi); % Finding frequency

% Eigenvector Normalization -- Unit Largest Entry:

new_U = [[0;0;U(:,8)]/max(abs(U(:,8))),[0;0;U(:,7)]/max(abs(U(:,7))),...
    [0;0;U(:,6)]/max(abs(U(:,6)))];
% ^^^^^^^^^^^^^^
% Needed to reverse the order of new_U to produce plots in the order of
% ascending frequencies from top (smallest) to bottom (largest)

%modes = U(:,6:8); % Modes

%% Plot Eigenvectors
for i = 1:3
    if i == 1
        col = 'b';
    elseif i ==2
        col = 'r';
    else
        col = 'g';
    end
    
    ploteigenvector(L,new_U(:,i),4,10,50,col,i)
end
%ploteigenvector(L,new_U(:,2),4,10,1,'r')
%ploteigenvector(L,new_U(:,3),4,10,1,'g')

%% Two-element model

%%Variables (English units)
L = 12; %in - Length of back half of plane
L_E = 4.5; %in - Elevator span
L_R = 5; %in - Rudder span
w = 1; %in - Width
h = 1/8; %in - Thickness of fuselage
h_E = 1/4; %in - Thickness of elevator
h_R = 0.04; %in - Thickness of rudder
E = 10175000; %psi - Elastic modulus
rho = 0.0002505; %lb-sec^2/in^4 - density
M_T = 1.131*rho; %[rho]*in^3 - Mass of tail
S_T = 0.5655*rho; %[rho]*in^4 - First mass-moment of tail
I_T = 23.124*rho; %[rho]*in^5 - Second mass-moment of tail
A = w*h;
I_zz = w*h^3/12;
cm2 = rho*A*L/100800;
ck2 = 4*E*I_zz/L^3;
cm4 = rho*A*L/806400;
ck4 = 8*E*I_zz/L^3;


M2 = cm2*[38544,0,5928,-642*L;0,344*L^2,642*L,-73*L^2;...
    5928,642*L,19272,-1458*L;-642*L,-73*L^2,-1458*L,172*L^2]+ ...
    [0,0,0,0;0,0,0,0;0,0,M_T,S_T;0,0,S_T,I_T];
K2 = ck2*[48,0,-24,6*L;0,4*L^2,-6*L,L^2;-24,...
    -6*L,24,-6*L;6*L,L^2,-6*L,2*L^2];

[eigvec_raw,omegasq,eigvec2_raw2] = eig(K2,M2);
omega2 = sqrt(diag(omegasq));
freq2 = omega2./(2*pi);

freq2 = freq2(1:3);

%validate
% syms w2
% omega22 = double(solve(det(K2-w2^2*M2) == 0));

%Change the eigenvectors (add the two zeros)
eigvec2 = [zeros(2,3);eigvec_raw(:,1:3)];

for i = 1:3
    if i == 1
        col = 'y';
    elseif i ==2
        col = 'm';
    else
        col = 'c';
    end
    
    ploteigenvector(L,eigvec2(:,i),2,40,1,col,i)
    if i == 1
        figure(1)
        legend('Mode 1, 4 element','Mode 1, 2 element')
    end
    if i == 2
        figure(1)
        legend('Mode 2, 4 element','Mode 2, 2 element')
    end
    if i == 3
        figure(1)
        legend('Mode 3, 4 element','Mode 3, 2 element')
    end
end