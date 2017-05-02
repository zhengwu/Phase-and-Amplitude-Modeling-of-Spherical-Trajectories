function R=AxisAndAngleRotation(u,theta)

ucross=[0 -u(3) u(2);
        u(3) 0 -u(1);
        -u(2) u(1) 0];
    

utensor=u*u';

R=eye(3)*cos(theta)+sin(theta)*ucross+(1-cos(theta))*utensor;

