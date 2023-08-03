% Adaptive Control - Simulation 2
% Masoud Pourghavam
% Student Number: 810601044
% Question 1-1

%% --------------------------------------------- %%

% Clear workspace, close figures, and clear command window
clear all;
close all;
clc;

%% System setup
format long;
s = tf('s');
g = 2*(s + 4)*(s + 0.58)/((s + 1.6)*(s + 3)*(s - 0.3));

%% Root locus plot
figure;
rlocus(g);
title('Root Locus');

%% System discretization and step response
l = feedback(g, 1);
w_b = bandwidth(l);
Ts = (2*pi)/(10*w_b);
Gz = c2d(g, Ts, 'zoh');

figure;
step(Gz);
title('Step Response');

%% Bode plot and pole-zero analysis
[numD, denD] = tfdata(Gz);
numD = cell2mat(numD);
denD = cell2mat(denD);

figure;
bode(g);
hold on;
bode(Gz);
title('Bode Plot');

poles = pole(Gz);
zeros = zero(Gz);

%% Desired pole location
Mp = 0.1;
zeta = ((log(Mp)^2)/(pi^2 + log(Mp)^2))^0.5;
ts = 3;
sigma = 4/ts;
wn = sigma/zeta;
s1 = -zeta*wn + 1i*(wn*(1-zeta^2)^0.5);
s2 = -zeta*wn - 1i*(wn*(1-zeta^2)^0.5);
s3 = -8*zeta*wn;
Pz1 = exp(s1*Ts);
Pz2 = exp(s2*Ts);
Pz3 = exp(s3*Ts);
z = tf('z');
Amm = (s - s1)*(s - s2)*(s - s3);
Am = (z - Pz1)*(z - Pz2)*(z - Pz3);

%% Calculations without zero cancellation
degA = length(denD) - 1;
A = 0;
for j = 1:length(denD)
    A = denD(j)*z^(length(denD)-j) + A;
end

degB = 0;
for j = 1:length(numD)
    if numD(j) ~= 0
        degB = degB + 1;
    end
end
degB = degB - 1;

B = 0;
for j = 1:length(numD)
    B = numD(j)*z^(length(numD)-j) + B;
end

%% Additional calculations
degB_plus = 0;
degB_minus = degB - degB_plus;
B_minus = B;
B_plus = 1;
degBm = degB;
deg_Bprim_m = degBm - degB_minus;
alpha = evalfr(Am, 1)/evalfr(B_minus, 1);
Bprim_m = alpha;
Bm = B_minus * Bprim_m;
deg_A0 = degA - degB_plus - 1;
A0 = z^(deg_A0);
T = A0 * Bprim_m;
degAm = degA;
degR = degAm - 1;
degR_prim = degR - degB_plus;

% Coefficient calculations
coefA = cell2mat(tfdata(A));
coefB_minus = cell2mat(tfdata(B_minus));
coefAm = cell2mat(tfdata(Am));
left = [1 0 coefB_minus(1) 0 0; 
        coefA(2) 1 coefB_minus(2) coefB_minus(1) 0;
        coefA(3) coefA(2) coefB_minus(3) coefB_minus(2) coefB_minus(1);
        coefA(4) coefA(3) 0 coefB_minus(3) coefB_minus(2);
        0 coefA(4) 0 0 coefB_minus(3)];
right = [coefAm(2) - coefA(2);
         coefAm(3) - coefA(3);
         coefAm(4) - coefA(4);
         0;
         0];
coefrs = (left^-1) * right;

R_prim = z^2 + coefrs(1) * z + coefrs(2);
S = coefrs(3) * z^2 + coefrs(4) * z + coefrs(5);
R = R_prim * B_plus;
G_cl = (B * T) / (A * R + B * S);

% Step response of closed-loop system
figure;
gg = minreal(G_cl);
step(gg);
title('Closed-Loop Step Response');
