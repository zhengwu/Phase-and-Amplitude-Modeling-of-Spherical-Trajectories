function q = path_to_q(p)

T = size(p,2);

for i=1:T-1    
    q(:,i)=(T-1)*InverseExp_Sphere(p(:,i),p(:,i+1));
    if norm(q(:,i))>0.0001
        q(:,i)=q(:,i)/sqrt(norm(q(:,i),'fro'));   
    else
        q(:,i)=[0;0;0];
    end
end
q(:,T) = ParallelTransport(p(:,T-1),q(:,T-1),p(:,T));


for i=2:T
    q(:,i)=BackwardParallelTranslation(p(:,1:i),q(:,i));
end

