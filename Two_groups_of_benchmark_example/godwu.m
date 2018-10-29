clear all;

[lambda,phi1,phi2] = PetscBinaryRead('wu_testtwogroups', 'indices','int32','precision','float64');
n=20;
m=24;
for i=1:m
    for j=1:n
        k=(i-1)*n+j;
        a(i,j) = phi1(k);
    end
end
for i=1:m
    for j=1:n
        k=(i-1)*n+j;
        b(i,j) = phi2(k);
    end
end

X=1:1:20
Y=1:1:24
figure(1)
surf(X,Y,a)
figure(2)
surf(X,Y,b)



