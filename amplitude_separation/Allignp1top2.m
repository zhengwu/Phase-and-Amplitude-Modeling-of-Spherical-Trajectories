function [pn1,indx,gam]=Allignp1top2(p1,p2)

T=size(p1,2);
N=21;

q1=path_to_q(p1);
q2=path_to_q(p2);


u1=p1(:,1);
u2=p2(:,1);

v1=InverseExp_Sphere(u1,p1(:,2));
v1=v1/norm(v1);
v2=InverseExp_Sphere(u2,p2(:,2));
v2=v2/norm(v2);

w1=cross(u1,v1);
w1=w1/norm(w1);
w2=cross(u2,v2);
w2=w2/norm(w2);

f1=[u1,v1,w1];
f2=[u2,v2,w2];

if ( norm(u1-u2)<10^-5 )
    % The starting points are the same
    %compute the cost
    len=0;
    t = linspace(0,1,T);
    [G,T1] = DynamicProgrammingQ2(q2,t,q1,t,t,t,0);
    gam = interp1(T1,G,t);
    q1n = Group_Action_by_Gamma_q(q1,gam);
    %     q2n = Group_Action_by_Gamma(q2,gamI);
    dmin=len^2+trapz(linspace(0,1,T),  sum( (q1n-q2).*(q1n-q2) )  );
    indx = 0;
    pn1 = q_to_path(q1n,u1);
else  
    
    M=120;
    theta=linspace(0,2*pi,M);
    for k=1:M
        R=AxisAndAngleRotation(u2,theta(k));
        v2R=R*v2;
        w2R=cross(u2,v2R);w2R=w2R/norm(w2R);
        f2=[u2,v2R,w2R];
        A=real(logm(f2*f1'));
        %compute \beta for each \theta
        for j=1:N
            B(:,j)=expm((j-1)*A/(N-1))*u1;
        end
        %parallel transpot of q1 along \beta to tangent space at u2
        for i=1:T
            q1par(:,i)=ForwardParallelTranslation(B,q1(:,i));
        end
        
        %compute the cost
        len=LengthOfTrajectory(B);
        t = linspace(0,1,T);
        [G,T1] = DynamicProgrammingQ2(q2,t,q1par,t,t,t,0);
        gam = interp1(T1,G,t);
        q1parn = Group_Action_by_Gamma_q(q1par,gam);
        TotalCost(k)=len^2+trapz(linspace(0,1,T),  sum( (q1parn-q2).*(q1parn-q2) )  );
        Part1(k) = len^2;
        Part2(k) = trapz(linspace(0,1,T),  sum( (q1parn-q2).*(q1parn-q2) )  );
    end
   
    [dmin,indx]=min(TotalCost);
    R=AxisAndAngleRotation(u2,theta(indx));
    v2R=R*v2;
    w2R=cross(u2,v2R);w2R=w2R/norm(w2R);
    f2=[u2,v2R,w2R];
    A=logm(f2*f1');
    %compute \beta for each \theta
    for j=1:N
        B(:,j)=expm((j-1)*A/(N-1))*u1;
    end
    
    %parallel transpot of q1 along \beta to tangent space at u2
    for i=1:T
        q1par(:,i)=ForwardParallelTranslation(B,q1(:,i));
    end
    
    [G,T1] = DynamicProgrammingQ2(q2,t,q1par,t,t,t,0);
    gam = interp1(T1,G,t);
    q1parn = Group_Action_by_Gamma_q(q1par,gam);
    
    for i=1:T
        q1n(:,i) = BackwardParallelTranslation(B,q1parn(:,i));
    end;
   
    pn1 = q_to_path(q1n,u1);

end;





