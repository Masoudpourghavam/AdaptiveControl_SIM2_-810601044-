function D = pluse(A,B)

if size(A,1) > size(A,2)
    A = A';
end

if size(B,1) > size(B,2)
    B = B';
end

na = length(A);
nb = length(B);

if na-nb > 0 
    A1 = A;
    B1 = [zeros(1,na-nb) B];
else
    A1 = [zeros(1,nb-na) A];
    B1 = B;
end


D1 = A1+B1;

i = 1;
while abs(D1(i)) <1e-7
    i = i+1;
end

D = D1(i:end);