function [R,S,E] = dioph(A1,B1,Ac1)

% clc
% clear 
% close

j = 1;
while A1(j) == 0
    j = j+1;
end
A = A1(j:end);
n = length(A)-1;

i = 1;
while B1(i) == 0
    i = i+1;
end
B2 = B1(i:end);

Ac = [zeros(1,2*n-length(Ac1)) Ac1];
B = [zeros(1,length(A)-length(B2)) B2];

% syms a1 a2 a3
% syms b0 b1 b2 b3
% syms al0 al1 al2 al3 al4 al5
% A = [1 a1 a2 a3];
% B = [b0 b1 b2 b3];
% n = length(A)-1;
% Ac = [al0 al1 al2 al3 al4 al5];

for j = 1:n
    E1(j,:) = [A(n-j+2:end) zeros(1,n-j)];
    E2(j,:) = [B(n-j+2:end) zeros(1,n-j)];
    E3(j,:) = [zeros(1,j-1) A(1:n-j+1)]; 
    E4(j,:) = [zeros(1,j-1) B(1:n-j+1)];
end

E = [E1 E2 ; E3 E4];

R_S = E\flipud(transpose(Ac));

R1 = flipud(R_S(1:n));
S1 = flipud(R_S(n+1:end));

k = 1;
while abs(R1(k)) <1e-3
    k = k+1;
end
R = R1(k:end);

d = 1;
while S1(d) == 0
    d = d+1;
end
S = S1(d:end);




