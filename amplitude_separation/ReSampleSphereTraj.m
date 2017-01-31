
function Xn = ReSampleSphereTraj(X,N)

[d,T] = size(X);
tau = 0:1/(T-1):1;
s = 0:1/(N-1):1;
Xn(:,1) = X(:,1)/norm(X(:,1));
Xn(:,N) = X(:,end)/norm(X(:,end));

for i=2:length(s)-1
    k=find(s(i)<=tau(:));
    ind1=k(1)-1;ind2=k(1);
    if(ind1==0)
        ind1=1;ind2=2;
    end
    w1=(s(i)-tau(ind1))/(tau(ind2)-tau(ind1));
    w2=(tau(ind2)-s(i))/(tau(ind2)-tau(ind1));
    p1 = X(:,ind1);
    p2 = X(:,ind2);
    
    theta = acos(sum(p1.*p2));
    Xn(:,i) = (1/sin(theta))* (sin(w2*theta)*p1+sin(w1*theta)*p2);
    Xn(:,i) = Xn(:,i)/norm(Xn(:,i)); 
end;