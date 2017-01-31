% inverse exponential on sphere

function v=InverseExp_Sphere(x,y)

theta=acos(sum(x.*y));
if theta <0.00001
    v=[0;0;0];
else
    v=theta/sin(theta)*( y-cos(theta)*x );
end

return;
