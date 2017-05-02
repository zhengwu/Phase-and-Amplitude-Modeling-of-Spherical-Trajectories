function qp=BackwardParallelTranslation(p,q)


T=size(p,2);

qp=q;

for i=T:-1:2
    qp=ParallelTransport(p(:,i),qp,p(:,i-1));
end

