
% Adaptive Control - Simulation 2
% Masoud Pourghavam
% Student Number: 810601044
% Question 3

%% --------------------------------------------- %%

clear all;
close all;
clc;

%%
num = conv([0.3 2],[0.95,-2]);
den = conv(conv([0.9 3],[0.82 -2]),[2 1]);
H = tf(num,den);
BW = bandwidth(H);
Ts = 0.1*pi/(BW);

%% " Indirect " STR " Without " Zero Cancellation

N = 5000;
landa = 1;
% ************** u_c **************
T = N/4;    % period
U = [ones(1,T/2) -ones(1,T/2)];
for k = 1:round(log(N)/log(2));
    U = [U U];
end
u_c = U(1:N)';
%  u_c = sin(.5*[1:N])';
%
%  u_c = (2/N)*[1:N]';

% % *******************************************
% % ************** Sys Definition *************

h = c2d(H,Ts,'zoh');
A = cell2mat(get(h,'den'));
B1 = cell2mat(get(h,'num'));
B = B1(2:end);

d0 = length(A)-length(B);
A_m = poly([0.5 0.6 -0.1]);
B_m = (sum(A_m)/sum(B))*B;

sys_m = idpoly([A_m],[B_m],[],[],[],0,Ts);
y_m = idsim([u_c],sys_m);

A_o = conv([1 .1],[1 -0.2])%[1 .5 0];
A_c = conv(A_o,A_m);

n = length(A)-1;
m = length(B)-1;

teta(:,3) = [1 2 1 1 1 -1]';
% teta(:,6) = randn(1,6);

p(:,1:6) = 1e6*eye(6);
p(1,1)=1e7;p(2,2)=1e7;p(3,3)=1e7;p(4,4)=1e5;p(5,5)=1e5;p(6,6)=1e5;

u(1:3) = zeros(1,3);
y(1:3) = zeros(1,3);

for t = 4:N
    
    y(t) = -A(2)*y(t-1) - A(3)*y(t-2) - A(4)*y(t-3) + B(1)*u(t-1) + B(2)*u(t-2) + B(3)*u(t-3);
    
    Fi = [-y(t-1) -y(t-2) -y(t-3) u(t-1) u(t-2) u(t-3)]';
    x = ( landa + Fi' * p(:,1:6) * Fi );
    epsil(t) = y(t) - Fi'*teta(:,t-1);
    p(:,7:12) = (p(:,1:6) - ( p(:,1:6) * Fi * Fi' * p(:,1:6) ) / x)/landa;
    K = p(:,1:6) * Fi / x;
    teta(:,t) = teta(:,t-1) + K*epsil(t);
    p(:,1:6) = p(:,7:12);
    
    A_hat = [1 teta(1:3,t)'];
    B_hat = teta(4:6,t)';
    
    [R,S] = dioph(A_hat,B_hat,A_c);
    T=(sum(A_m)/sum(B_hat))*A_o;
    T11=(sum(A_m)/sum(B))*A_o;
    for j=1:3
        TT(t-3,j) = T(j);
        TT11(t-3,j)=T11(j);
    end
    u(t) = -R(2)*u(t-1) -R(3)*u(t-2) +T(1)*u_c(t) +T(2)*u_c(t-1) +T(3)*u_c(t-2) -S(1)*y(t) -S(2)*y(t-1) -S(3)*y(t-2);
end
ff=TT-TT11;
figure(1)
plot(ff(1:100,1),'b')
hold on
plot(ff(1:100,2),'k')
hold on 
plot(ff(1:100,3),'--r')
title('T Parameters Estimation Error')
legend('ET_0','ET_1','ET_2')


A_close = pluse(conv(A_hat,R),conv(B_hat,S));
B_y = conv(B_hat,T);
B_u = conv(A_hat,T);
figure(2)
plot(u_c,'--r','LineWidth',2)
hold on
%plot(y_m,'g')
plot(y,'b','LineWidth',1.5)
axis([0 N/2 -1.5 1.5])

title('Output - Indi. STR without zero cancellation for N.M.P. sys.')
legend('u_c','y')


figure(3)
plot(u,'LineWidth',2)
axis([0 N/2 -100 100])
ylim([-5,5])
title('Control Signal - Indi. STR without zero cancellation for N.M.P Sys.')


figure(4)
subplot(6,1,1), plot(A(2)*ones(1,N),'r--','LineWidth',1.5), hold on, plot(teta(1,:),'LineWidth',1.4),xlim([0 N/20]),ylabel('a_1')
subplot(6,1,2), plot(A(3)*ones(1,N),'r--','LineWidth',1.5), hold on, plot(teta(2,:),'LineWidth',1.4),xlim([0 N/20]),ylabel('a_2')
subplot(6,1,3), plot(A(4)*ones(1,N),'r--','LineWidth',1.5), hold on, plot(teta(3,:),'LineWidth',1.4),xlim([0 N/20]),ylabel('a_3')
subplot(6,1,4), plot(B(1)*ones(1,N),'r--','LineWidth',1.5), hold on, plot(teta(4,:),'LineWidth',1.4),xlim([0 N/20]),ylabel('b_1')
subplot(6,1,5), plot(B(2)*ones(1,N),'r--','LineWidth',1.5), hold on, plot(teta(5,:),'LineWidth',1.4),xlim([0 N/20]),ylabel('b_2')
subplot(6,1,6), plot(B(3)*ones(1,N),'r--','LineWidth',1.5), hold on, plot(teta(6,:),'LineWidth',1.4),xlim([0 N/20]),ylabel('b_3')


