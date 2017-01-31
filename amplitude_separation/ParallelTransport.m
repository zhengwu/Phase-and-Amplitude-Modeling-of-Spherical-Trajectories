function vtilde=ParallelTransport(x,v,y)

vtilde=v-2*sum(v.*y)*(x+y)/(norm(x+y,'fro'))^2;