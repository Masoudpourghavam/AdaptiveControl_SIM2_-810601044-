% Adaptive Control - Simulation 2
% Masoud Pourghavam
% Student Number: 810601044
% Question 2-4 without color noise

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


%
n=6;
a2=denD(2);
a1=denD(3);
a0=denD(4);
b2=numD(2);
b1=numD(3);
b0=numD(4);
teta0=[a2;a1;a0;b2;b1;b0];

N = 100;
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
 
y_out=zeros(N,1);
%initial condition
phi_out=zeros(1,n);
k=1;
sse=0;
I=0;
NN=6; % number of controller parameter
teta_hat=[0.1;0.1;0.1;0;0;0];
p=1000*eye(NN);
phi_t=zeros(1,NN);

r_hat2=zeros(N,1);
r_hat1=zeros(N,1);
r_hat0=zeros(N,1);
s_hat2=zeros(N,1);
s_hat1=zeros(N,1);
s_hat0=zeros(N,1);


u_control=zeros(N,1);

RR=zeros(1,2);
SS=zeros(1,3);

Am=(z-Pz1)*(z-Pz2)*(z-Pz3);
Am_coef=cell2mat(tfdata(Am));
degAm=3;
degBm=2;
degB_plus=degBm;
degB_minus=degBm-degB_plus;
d0=1;

bm0=evalfr(Am,1);
Bm=bm0*z^(degBm);
% deg_Bprim_m=degBm-degB_minus;
% Bprim_m=(bm0/B_minus)*z^(deg_Bprim_m);
deg_A0=degAm-degB_plus-1;
A0=z^(d0-1);
A0_star=A0*(z^(1-d0));
T=A0_star*evalfr(Am,1);
T_coef=cell2mat(tfdata(T));
t2_hat=T_coef(1);

yf=zeros(N,1);
uf=zeros(N,1);
yfs=zeros(1,3);
ufs=zeros(1,3);
%
for k=2:N-1
    
    for ii=1:3
        if k-ii<=0
            ud=0;
            yd=0;
        else
            ud=-uf(k-ii,1);
            yd=-yf(k-ii,1);
        end
        ufs(1,ii)=ud;
        yfs(1,ii)=yd;
    end
    uf(k,1)=u_control(k,1)+ufs*[Am_coef(2);Am_coef(3);Am_coef(4)];
    yf(k,1)=y_out(k,1)+yfs*[Am_coef(2);Am_coef(3);Am_coef(4)];
   
    for m=0:2
        if k-m-d0<=0
            y_dummy=0;
            u_dummy=0;
        else 
            y_dummy=yf(k-m-d0,1);
            u_dummy=uf(k-m-d0,1);
        end
        phi_t(1,m+1)=u_dummy;
        phi_t(1,m+4)=y_dummy;
    end
    p=p-((p*(phi_t')*phi_t*p)/(1+phi_t*p*(phi_t')));
    gain=p*(phi_t');
    
    for jj=0:2
        if k-jj<=0
            y_dummyy=0;
            u_dummyy=0;
        else 
            y_dummyy=-y_out(k-jj,1);
            u_dummyy=u_control(k-jj,1);
        end
        phi_out(1,jj+1)=y_dummyy;
        phi_out(1,jj+4)=u_dummyy;
    end
    
    y_out(k+1,1)=phi_out*teta0;
    
    teta_hat=teta_hat+gain*(y_out(k,1)-(phi_t*teta_hat));
    
    r_hat2(k,1)=teta_hat(1);
    r_hat1(k,1)=teta_hat(2);
    r_hat0(k,1)=teta_hat(3);   
    s_hat2(k,1)=teta_hat(4);
    s_hat1(k,1)=teta_hat(5);
    s_hat0(k,1)=teta_hat(6);
    
    for ii=0:1
        if k-ii<=0
            uu=0;
        else
            uu=-u_control(k-ii,1);
        end
        RR(1,ii+1)=uu;
    end
    for ii=0:2
        if k-ii<=0
            yy=0;
        else 
            yy=y_out(k-ii+1,1);
        end
        SS(1,ii+1)=yy;    
    end
    
    u_control(k+1,1)=(RR*[r_hat1(k,1);r_hat0(k,1)]+t2_hat*uc_in(k,1)-SS*[s_hat2(k,1);s_hat1(k,1);s_hat0(k,1)])/r_hat2(k,1);
   
    I=(abs(y_out(k,1)-(phi_t*teta_hat)));
    sse=I;
end
%%
teta0
teta_hat
%%
figure 
plot(sample_number,uc_in, "black", 'LineWidth', 2)
hold on 
plot(sample_number,y_out, "green", 'LineWidth', 2)
ylim([-2 2])
xlabel('Iter')
legend('Uc','y')

figure
plot(sample_number,u_control, "black", 'LineWidth', 2)
ylim([-50 50])
xlabel('Iter')
ylabel('u control')
%%
figure
subplot(3,1,1) 
plot(sample_number,r_hat2, "black", 'LineWidth', 2)
xlabel('Iter')
ylabel('rhat2')
subplot(3,1,2)
plot(sample_number,r_hat1, "green", 'LineWidth', 2)
xlabel('Iter')
ylabel('rhat1')
subplot(3,1,3)
plot(sample_number,r_hat0, "red", 'LineWidth', 2)
xlabel('Iter')
ylabel('rhat0')
%%
figure
subplot(3,1,1)
plot(sample_number,s_hat2, "black", 'LineWidth', 2)
xlabel('Iter')
ylabel('shat2')

subplot(3,1,2)
plot(sample_number,s_hat1, "green", 'LineWidth', 2)
xlabel('Iter')
ylabel('shat1')

subplot(3,1,3)
plot(sample_number,s_hat0, "red", 'LineWidth', 2)
xlabel('Iter')
ylabel('shat0')
