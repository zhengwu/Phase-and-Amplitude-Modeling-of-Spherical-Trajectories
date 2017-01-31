function [mu,psi,vec] = SqrtMean(gam)
 
[n,T] = size(gam);

psi = zeros(n,T-1);
for i=1:n
    psi(i,:) = sqrt(diff(gam(i,:))*T);
end

%Find direction
mu = psi(1,:);
t = 1;

clear vec;
for iter = 1:5
    for i=1:n
        v = psi(i,:) - mu;
        len = acos(sum(mu.*psi(i,:))/T);
        if len > 0.05
            vec(i,:) = (len/sin(len))*(psi(i,:) - cos(len)*mu);
        else
            vec(i,:) = zeros(1,T-1);
        end
    end
    vm = mean(vec);
    lvm(iter) = sqrt(sum(vm.*vm)/T);
    %[iter lvm(iter)]
    mu = cos(t*lvm(iter))*mu + (sin(t*lvm(iter))/lvm(iter))*vm;
end

for i=1:n
    phi(i,:) = cumsum([0 psi(i,:).*psi(i,:)/T]);
end
