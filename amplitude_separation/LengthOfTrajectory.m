function [ len ] = LengthOfTrajectory( p )

[n,N]=size(p);

for i=1:n
    v(i,:)=gradient(p(i,:),1/(N-1));
end

len=sqrt(trapz(linspace(0,1,N),sum(v.*v)));

end

