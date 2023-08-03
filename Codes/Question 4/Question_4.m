% Adaptive Control - Simulation 2
% Masoud Pourghavam
% Student Number: 810601044
% Question 4

%% --------------------------------------------- %%

clear all;
close all;
clc;

%%

s=tf('s');
g=4/((s+0.75)*(s+3));
n=2;
degA=2;
degB=0;
B=4;
A=((s+0.4)*(s+3.8));
Mp=0.1;
zeta=((log(Mp)^2)/(pi^2+log(Mp)^2))^0.5;
ts=3;
sigma=4/ts;
wn=sigma/zeta;
Bm=wn^2;
Am=s^2+(2*zeta*wn)*s+wn^2;

% without cancelation 
B_plus=1;
B_minus=B;
B_prim_m=Bm/B_minus;
degA0=degA-1;
a0_const=20;
A0=s+a0_const;
Hf=1/Am;
est_param=3;
p=1000*eye(est_param);
teta_hat=ones(est_param,1);
t0=B_prim_m;
alpha=0.3;
