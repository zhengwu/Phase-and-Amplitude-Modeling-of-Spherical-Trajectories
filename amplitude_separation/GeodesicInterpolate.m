function gam = GeodesicInterpolate(t0,t1,p0,p1)


d = acos(sum(p0.*p1));
if d < 0
    d= d + 2*pi;
end

if d==0
    for t=t0:1:t1
        i = t - t0;
        gam(:,i+1) = p0;
    end    
else
    for t=t0:1:t1
        i = t - t0;
        tau = (t - t0)/(t1-t0);
        gam(:,i+1) = (sin((1-tau)*d)*p0 + sin(tau*d)*p1)/sin(d);
    end
end


