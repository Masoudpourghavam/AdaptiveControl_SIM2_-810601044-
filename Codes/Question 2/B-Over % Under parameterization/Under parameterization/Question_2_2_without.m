% Adaptive Control - Simulation 2
% Masoud Pourghavam
% Student Number: 810601044
% Question 2-2 without color noise with under parameterization

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
n=4;  % number of system parameter
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

u_control=zeros(N,1);

Am=(z-Pz1)*(z-Pz2)*(z-Pz3);
Am_coef=cell2mat(tfdata(Am));
degAm=3;
degBm=2;
degB_plus=0;
degB_minus=degBm-degB_plus;
deg_Bprim_m=degBm-degB_minus;
B_plus=1;
B=0;
for j=1:length(numD)
   B=numD(j)*z^(length(numD)-j)+B; 
end
B_minus=B;
alpha=evalfr(Am,1)/evalfr(B_minus,1);
Bprim_m=alpha*z^(deg_Bprim_m);
Bm=B_minus*Bprim_m;
Bm_coef=cell2mat(tfdata(Bm));
d0=1;

deg_A0=degAm-degB_plus-1;
A0=z^(deg_A0);
A0_star=A0*(z^(-deg_A0));
T=A0*Bprim_m;
T_coef=cell2mat(tfdata(T));
t2_hat=T_coef(1);
degR=degAm-1;
degS=degAm-1;
degR_italic=degR+degB_minus;
degS_italic=degS+degB_minus;
Rhat_italic=zeros(N,degR_italic+1);
Shat_italic=zeros(N,degS_italic+1);

NN=degR_italic+degS_italic+2; % number of controller parameter
teta_hat=ones(NN,1)+rand(NN,1);
p=1000*eye(NN);
phi_t=zeros(1,NN);

r_hat=zeros(N,degR+1);
s_hat=zeros(N,degS+1);


yf=zeros(N,1);
uf=zeros(N,1);

yfs=zeros(1,3);
ufs=zeros(1,3);
syms sss
roots_Ritalic=zeros(1,degR_italic);
roots_Sitalic=zeros(1,degS_italic);
diffrence_matrix=zeros(degR_italic,degS_italic);
eps=0.1;
vv=0;
ee=0;
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
   
    for m=0:degR_italic
        if k-m-d0<=0
            y_dummy=0;
            u_dummy=0;
        else 
            y_dummy=yf(k-m-d0,1);
            u_dummy=uf(k-m-d0,1);
        end
        phi_t(1,m+1)=u_dummy;
        phi_t(1,m+degR_italic+2)=y_dummy;
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
    
    Rhat_italic(k,:)=teta_hat(1:degR_italic+1);
    Shat_italic(k,:)=teta_hat(degR_italic+2:end);

   
	R_italic=0;
	S_italic=0;
	for kk=0:degR_italic
        R_italic=R_italic+Rhat_italic(k,degR_italic+1-kk)*sss^(kk); 
	end
	for kk=0:degS_italic
        S_italic=S_italic+Shat_italic(k,degS_italic+1-kk)*sss^(kk);
	end
	R_new=R_italic;
	S_new=S_italic;
	roots_Ritalic=double(vpa(root(R_italic)));
	roots_Sitalic=double(vpa(root(S_italic))); 
	for ii=1:degR_italic
        for jj=1:degS_italic
          diffrence_matrix(ii,jj)=abs(roots_Ritalic(ii)-roots_Sitalic(jj));
        end
	end
   
	for ii=1:degR_italic
        for jj=1:degS_italic
            if diffrence_matrix(ii,jj)<=eps
           %   if ((imag(roots_Ritalic(ii))==0 && imag(roots_Sitalic(jj)==0)) || (imag(roots_Ritalic(ii))~=0 && imag(roots_Sitalic(jj)~=0)))
                 sd=sss-roots_Sitalic(jj);
%                 
                 rd=sss-roots_Ritalic(ii);
%                
                 
                 [qq,rr]=quorem(R_new,rd);
                 R_new=qq;
                 
                 [qqs,rrs]=quorem(S_new,sd);
                 S_new=qqs;
                
          %    end
            end
        end
	end
    
	aa=double(vpa(coeffs(R_new,'All')));
 
	if length(aa)==degR+1
        vv=vv+1;
        r_hat(vv,:)=aa;
	end
	bb=double(vpa(coeffs(S_new,'All')));
	if length(bb)==degS+1
        ee=ee+1;
        s_hat(ee,:)=bb; 
	end
 
	RR=zeros(1,length(aa)-1);
	SS=zeros(1,length(bb));

	for ii=0:length(aa)-2
        if k-ii<=0
            uu=0;
        else
            uu=-u_control(k-ii,1);
        end
        RR(1,ii+1)=uu;
	end
	for ii=0:length(bb)-1
        if k-ii<=0
            yy=0;
        else 
            yy=y_out(k-ii+1,1);
        end
        SS(1,ii+1)=yy;    
	end 
	u_control(k+1,1)=(RR*(aa(2:end)')+t2_hat*uc_in(k,1)-SS*(bb'))/aa(1);
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
legend ('u control')
%%
rhat2=r_hat(1:vv,1);
rhat1=r_hat(1:vv,2);
rhat0=r_hat(1:vv,3);
shat2=s_hat(1:ee,1);
shat1=s_hat(1:ee,2);
shat0=s_hat(1:ee,3);
%%
figure
subplot(3,1,1)
plot(1:vv,rhat2, "black", 'LineWidth', 2)
xlabel('Iter')
ylabel('rhat2')
legend ('rhat2')
subplot(3,1,2)
plot(1:vv,rhat1, "green", 'LineWidth', 2)
xlabel('Iter')
ylabel('rhat1')
legend ('rhat1')
subplot(3,1,3) 
plot(1:vv,rhat0, "red", 'LineWidth', 2)
xlabel('Iter')
ylabel('rhat0')
legend ('rhat0')
%%
figure 
subplot(3,1,1)
plot(1:ee,shat2, "black", 'LineWidth', 2)
xlabel('Iter')
ylabel('shat2')
legend ('shat2')
subplot(3,1,2) 
plot(1:ee,shat1, "green", 'LineWidth', 2)
xlabel('Iter')
ylabel('shat1')
legend ('shat1')
subplot(3,1,3)
plot(1:ee,shat0, "red", 'LineWidth', 2)
xlabel('Iter')
ylabel('shat0')
legend ('shat0')

