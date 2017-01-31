function qp=ForwardParallelTranslation(p,q)
%p: path
%q: vector field needed to be translated.

T=size(p,2);

qp=q;

for i=1:T-1
    qp=ParallelTransport(p(:,i),qp,p(:,i+1));
end

