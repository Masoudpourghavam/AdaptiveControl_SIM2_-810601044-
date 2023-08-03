% Adaptive Control - Simulation 2
% Masoud Pourghavam
% Student Number: 810601044
% Question 2-1 With color noise with over parameterization

%% --------------------------------------------- %%

clear all;
close all;
clc;

%%

format long
s=tf('s');
g=2*(s+4)*(s+0.58)/((s+1.6)*(s+3)*(s-0.3));
l=feedback(g,1);
w_b=bandwidth(l);
Ts=(2*pi)/(10*w_b);
Gz=c2d(g,Ts,'zoh');
figure
step(Gz)
[numD,denD]=tfdata(Gz);
numD=cell2mat(numD);
denD=cell2mat(denD);
figure
bode(g)
hold on
bode(Gz)
pole=pole(Gz);
zero=zero(Gz);
%% desired pole location 
Mp=0.1;
zeta=((log(Mp)^2)/(pi^2+log(Mp)^2))^0.5;
ts=3;
sigma=4/ts;
wn=sigma/zeta;
s1=-zeta*wn+i*(wn*(1-zeta^2)^0.5);
s2=-zeta*wn-i*(wn*(1-zeta^2)^0.5);
s3=-8*zeta*wn;
Pz1=exp(s1*Ts);
Pz2=exp(s2*Ts);
Pz3=exp(s3*Ts);
z=tf('z');
Am=(z-Pz1)*(z-Pz2)*(z-Pz3);
%%
n=8;
a2=denD(2);
a1=denD(3);
a0=denD(4);
b2=numD(2);
b1=numD(3);
b0=numD(4);
teta0=[a2;a1;a0;b2;b1;b0];
N = 500;
sample_number=zeros(N,1);
for i=1:N
    sample_number(i,1)=i;
end

%square pulse refrence input
uc_in=zeros(N,1);
for i=1:N/4
    uc_in(i,1)=1;
end
for i=N/4:N/2
    uc_in(i,1)=-1;
end
for i=N/2:3*N/4
    uc_in(i,1)=1;
end
for i=3*N/4:N
    uc_in(i,1)=-1;
end
figure
plot(sample_number,uc_in)

% calculate output with noise 
y_out=zeros(N,1);
epsilon=zeros(N,1);
%Generate color noise 
stdv_noise=0.05;
noise=stdv_noise*rand(N,1);
noise=noise-mean(noise);
c0=0.5;
for i=2:N 
        noise(i)=c0*noise(i-1)+noise(i);
end
teta0=[teta0;c0];

%initial condition
NN=7;
teta_hat=[0.01;0.02;0.03;0.04;0.05;0.06;0.1];
p=1000*eye(NN);
phi_t=zeros(1,NN);
k=1;
sse=0;
I=0;

a_hat2=zeros(N,1);
a_hat1=zeros(N,1);
a_hat0=zeros(N,1);
b_hat2=zeros(N,1);
b_hat1=zeros(N,1);
b_hat0=zeros(N,1);
c_hat0=zeros(N,1);
r1=zeros(N,1);
r0=zeros(N,1);
t0=zeros(N,1);
s2=zeros(N,1);
s1=zeros(N,1);
s0=zeros(N,1);
u_control=zeros(N,1);
RR=zeros(1,2);
SS=zeros(1,3);
%%
for k=1:N
    for m=1:3
        if k-m<=0
            y_dummy=0;
            u_dummy=0;
            e_dummy=0;
        else 
            y_dummy=-y_out(k-m,1);
            u_dummy=u_control(k-m,1);
            e_dummy=noise(k-m,1);
        end
        phi_t(1,m)=y_dummy;
        phi_t(1,m+3)=u_dummy;
        phi_t(1,end)=e_dummy;
    end
    p=p-((p*(phi_t')*phi_t*p)/(1+phi_t*p*(phi_t')));
    gain=p*(phi_t');
    y_out(k,1)=phi_t*teta0+noise(k,1);
    epsilon(k,1)=y_out(k,1)-phi_t*teta_hat;
    teta_hat=teta_hat+gain*(epsilon(k,1));
    a_hat2(k,1)=teta_hat(1);
    a_hat1(k,1)=teta_hat(2);
    a_hat0(k,1)=teta_hat(3);   
    b_hat2(k,1)=teta_hat(4);
    b_hat1(k,1)=teta_hat(5);
    b_hat0(k,1)=teta_hat(6);
    c_hat0(k,1)=teta_hat(7);
    A_hat=z^3+a_hat2(k,1)*z^2+a_hat1(k,1)*z+a_hat0(k,1);
    B_hat=b_hat2(k,1)*z^2+b_hat1(k,1)*z+b_hat0(k,1);
    
    degA_hat=3;
    degB_hat=2;
    degB_plus=0;
    degB_minus=degB_hat-degB_plus;
    B_minus=B_hat;
    B_plus=1;
    degBm=degB_hat;
    alpha=evalfr(Am,1)/evalfr(B_minus,1);
    deg_Bprim_m=degBm-degB_minus;
    Bprim_m=alpha*z^(deg_Bprim_m);
    Bm=B_minus*Bprim_m;
    
    deg_A0=degA_hat-degB_plus-1;
    A0=z^(deg_A0);
    
    T=A0*Bprim_m;
    Tcoef=cell2mat(tfdata(T));
    t0(k,1)=Tcoef(1);
    
    
coefA_hat=cell2mat(tfdata(A_hat));
coefB_minus=cell2mat(tfdata(B_minus));
coefAm=cell2mat(tfdata(Am));
left=[1 0 coefB_minus(1) 0 0;coefA_hat(2) 1 coefB_minus(2) coefB_minus(1) 0;...
    coefA_hat(3) coefA_hat(2) coefB_minus(3) coefB_minus(2) coefB_minus(1);...
    coefA_hat(4) coefA_hat(3) 0 coefB_minus(3) coefB_minus(2);...
    0 coefA_hat(4) 0 0 coefB_minus(3)];
right=[coefAm(2)-coefA_hat(2);coefAm(3)-coefA_hat(3);coefAm(4)-coefA_hat(4);0;0];
coefrs=(left^-1)*right;
    
    
    R_prim=[1 coefrs(1) coefrs(2)];
    R=R_prim*B_plus;
    r1(k,1)=R(2);
    r0(k,1)=R(3);
    

    Scoef=[coefrs(3) coefrs(4) coefrs(5)];
    s2(k,1)=Scoef(1);
    s1(k,1)=Scoef(2);
    s0(k,1)=Scoef(3);
    
    for ii=1:2
        if k-ii<=0
            uu=0;
        else
            uu=-u_control(k-ii,1);
        end
        RR(1,ii)=uu;
    end
    for ii=0:2
        if k-ii<=0
            yy=0;
        else 
            yy=y_out(k-ii,1);
        end
        SS(1,ii+1)=yy;    
    end
    
    u_control(k,1)=RR*[r1(k,1);r0(k,1)]+t0(k,1)*uc_in(k,1)-SS*[s2(k,1);s1(k,1);s0(k,1)];
   
    I=(abs(y_out(k,1)-(phi_t*teta_hat)));
    sse=I;
end
%%
figure 
plot(sample_number,uc_in, "green", 'LineWidth', 2)
hold on 
plot(sample_number,y_out, "black", 'LineWidth', 2)
ylim([-2 2])
xlabel('Iter')
legend('Uc','y')
%%
figure
plot(sample_number,u_control, "black", 'LineWidth', 2)
ylim([-50 50])
xlabel('Iter')
ylabel('u control')
legend('u')
%%
figure
subplot(3,1,1) 
plot(sample_number,a_hat2, "black", 'LineWidth', 2)
xlabel('Iter')
ylabel('ahat2')
legend('ahat2')
subplot(3,1,2)
plot(sample_number,a_hat1, "green", 'LineWidth', 2)
xlabel('Iter')
ylabel('ahat1')
legend('ahat1')
subplot(3,1,3)
plot(sample_number,a_hat0, "red", 'LineWidth', 2)
xlabel('Iter')
ylabel('ahat0')
legend('ahat0')
%%
figure
subplot(3,1,1)
plot(sample_number,b_hat2, "black", 'LineWidth', 2)
xlabel('Iter')
ylabel('bhat2')
legend('bhat2')
subplot(3,1,2)
plot(sample_number,b_hat1, "green", 'LineWidth', 2)
xlabel('Iter')
ylabel('bhat1')
legend('bhat1')
subplot(3,1,3)
plot(sample_number,b_hat0, "red", 'LineWidth', 2)
xlabel('Iter')
ylabel('bhat0')
legend('bhat0')
%%
figure
plot(sample_number,c_hat0, "black", 'LineWidth', 2)
xlabel('Iter')
ylabel('chat0')
%%
figure
subplot(2,1,1)
plot(sample_number,r0, "black", 'LineWidth', 2)
xlabel('Iter')
ylabel('r0')
legend('r0')
subplot(2,1,2)
plot(sample_number,r1, "green", 'LineWidth', 2)
xlabel('Iter')
ylabel('r1')
legend('r1')
%%
figure
subplot(3,1,1)
plot(sample_number,s0, "black", 'LineWidth', 2)
xlabel('Iter')
ylabel('s0')
legend('s0')
subplot(3,1,2)
plot(sample_number,s1, "green", 'LineWidth', 2)
xlabel('Iter')
ylabel('s1')
legend('s1')
subplot(3,1,3)
plot(sample_number,s2, "red", 'LineWidth', 2)
xlabel('Iter')
ylabel('s2')
legend('s2')
%%
figure
plot(sample_number,t0, "black", 'LineWidth', 2)
xlabel('Iter')
ylabel('t0')
legend('t0')
%%
figure
Rc=R(1)*z^2+R(2)*z+R(3);
S=Scoef(1)*z^2+Scoef(2)*z+Scoef(3);
G_cl=(B_hat*T)/(A_hat*Rc+B_hat*S);
step(G_cl);