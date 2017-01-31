function gam = PiecewiseGeodesic(T,tau,P)


k = size(P,2);

gam(:,1) = P(:,1);
if tau(1) > 1
      gam(:,1:tau(1)) = repmat(P(:,1),1,tau(1));
end
for i=1:k-1
    if tau(i+1) - tau(i) > 1
        gam(:,tau(i):tau(i+1)) = GeodesicInterpolate(tau(i),tau(i+1),P(:,i),P(:,i+1));
    end
end
if tau(k) < T
    gam(:,tau(k):T) = GeodesicInterpolate(tau(k),T,P(:,k),P(:,k));
end

for i=1:k
    gam(:,tau(i))=P(:,i);
end