function qn = Group_Action_by_Gamma_q(q,gamma)

[n,T] = size(q);
gammadot = gradient(gamma,1/T);
qn = zeros(n,T);
for j = 1:n
qn(j,:)=spline(linspace(0,1,T) , q(j,:) ,gamma).*sqrt(gammadot);
end;
