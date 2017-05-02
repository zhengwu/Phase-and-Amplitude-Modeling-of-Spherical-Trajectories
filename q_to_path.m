% q_to_path

function p=q_to_path(q,p0)

T=size(q,2);
stp=1.0/(T-1);
p(:,1)=p0;
p(:,2) = Exp_Sphere(p0,stp*q(:,1)*norm(q(:,1),'fro'));
for i=2:T-1
    v = ForwardParallelTranslation(p(:,1:i),q(:,i)*norm(q(:,i),'fro'));    
    p(:,i+1)=Exp_Sphere(p(:,i),stp*v);
%     p(:,i+1)=Exp_Sphere(p(:,i),stp*q(:,i)*norm(q(:,i),'fro'));
end

% project to the sphere;
for i=1:T
    p(:,i) = p(:,i)/norm(p(:,i),'fro');
end;