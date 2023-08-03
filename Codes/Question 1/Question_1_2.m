% Adaptive Control - Simulation 2
% Masoud Pourghavam
% Student Number: 810601044
% Question 1-2

%% --------------------------------------------- %%

% Clear workspace, close figures, and clear command window
clear all;
close all;
clc;

%% System setup
format long;
s = tf('s');
g = 2*(s + 4)*(s + 0.58)/((s + 1.6)*(s + 3)*(s - 0.3));

%% Step response of discretized system
l = feedback(g, 1);
w_b = bandwidth(l);
Ts = (2*pi)/(10*w_b);
Gz = c2d(g, Ts, 'zoh');
figure();
step(Gz);
title('Step Response of Discretized System');

%% Bode plot and pole-zero analysis
[numD, denD] = tfdata(Gz);
numD = cell2mat(numD);
denD = cell2mat(denD);

figure();
bode(g);
hold on;
bode(Gz);
legend('Continuous', 'Discrete');
grid on;
title('Bode Plot');

poles = pole(Gz);
zeros_Gz = zero(Gz);

%% Desired pole location
Mp = 0.1;
zeta = ((log(Mp)^2)/(pi^2 + log(Mp)^2))^0.5;
ts = 3;
sigma = 4/ts;
wn = sigma/zeta;
s1 = -zeta*wn + 1i*(wn*(1-zeta^2)^0.5);
s2 = -zeta*wn - 1i*(wn*(1-zeta^2)^0.5);
s3 = -8*zeta*wn;

Gm = 1/((s - s1)*(s - s2)*(s - s3));
figure();
step(Gm);
title('Step Response of Gm');

Gmz = c2d(Gm, Ts, 'zoh');
figure();
step(Gmz);
title('Step Response of Discretized Gm');

%% Am calculation
Pz1 = exp(s1*Ts);
Pz2 = exp(s2*Ts);
Pz3 = exp(s3*Ts);
z = tf('z');
Am = (z - Pz1)*(z - Pz2)*(z - Pz3);

%% Zero cancellation
degA = length(denD) - 1;
A = 0;
for j = 1:length(denD)
    A = denD(j)*z^(length(denD) - j) + A;
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
    B = numD(j)*z^(length(numD) - j) + B;
end

degB_plus = degB;
degB_minus = degB - degB_plus;
for j = 1:length(numD)
    if numD(j) ~= 0
        beta = numD(j);
        break
    end
end

B_minus = beta;
B_plus = B / B_minus;
degBm = degB;
bm0 = evalfr(Am, 1);
Bm = bm0 * z^(degBm);
deg_Bprim_m = degBm - degB_minus;
Bprim_m = (bm0 / beta) * z^(deg_Bprim_m);
deg_A0 = degA - degB_plus - 1;
A0 = z^(deg_A0);
R_prim = 1;
R = R_prim * B_plus;
T = A0 * Bprim_m;
S = (Am - A) / B_minus;
G_cl = (B * T) / (A * R + B * S);

% Step response of closed-loop system
figure();
step(G_cl);
title('Closed-Loop Step Response');

gg = minreal(G_cl);
figure();
step(gg);
title('Minimized Closed-Loop Step Response');
