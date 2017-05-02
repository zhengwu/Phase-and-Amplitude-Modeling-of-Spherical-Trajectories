function pn = Group_Action_by_Gamma_p(p,gamma)

[n,T] = size(p);
% for j = 1:n
% pn(j,:)=spline(linspace(0,1,T) , p(j,:) ,gamma);
% end;
% for i=1:T
%     pn(:,i)=pn(:,i)/norm(pn(:,i));
% end

tau = 0:1/(T-1):1;
s = gamma;
pn(:,1) = p(:,1);
pn(:,T) = p(:,end);

for i=2:length(s)-1
    k=find(s(i)<=tau(:));
    ind1=k(1)-1;ind2=k(1);
    if(ind1==0)
        ind1=1;ind2=2;
    end
    w1=(s(i)-tau(ind1))/(tau(ind2)-tau(ind1));
    w2=(tau(ind2)-s(i))/(tau(ind2)-tau(ind1));
    p1 = p(:,ind1);
    p2 = p(:,ind2);
    
    theta = acos(sum(p1.*p2));
    if(theta<10^-3)
        pn(:,i) = p1;
    else
        pn(:,i) = (1/sin(theta))* (sin(w2*theta)*p1+sin(w1*theta)*p2);
        pn(:,i) = pn(:,i)/norm(pn(:,i));
    end
end;
