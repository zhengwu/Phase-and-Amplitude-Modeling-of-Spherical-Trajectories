function [epu1,epw1,dsq]=InverseExponentialMap(p1,p2)
%q2 - q1


T=size(p1,2);
N=21;
disppic = 0;

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
    
if ( norm(u1-u2)<0.005 )
    
    % The starting points are the same
    %compute the cost
    len=0;
    t = linspace(0,1,T);
    [G,T1] = DynamicProgrammingQ2(q2,t,q1,t,t,t,0);
    gam = interp1(T1,G,t);
    q1n = Group_Action_by_Gamma_q(q1,gam);
    dsq=len^2+trapz(linspace(0,1,T),  sum( (q1n-q2).*(q1n-q2) )  );
    epw1 = q2 - q1n;
    epu1 = zeros(3,1);
    
else    
    %inital theta
    M=60;
    theta=linspace(0,2*pi,M);
    oldgam = zeros(1,50);
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
        
        % with registration
        t = linspace(0,1,T);
        [G,T1] = DynamicProgrammingQ2(q2,t,q1par,t,t,t,0);
        gam = interp1(T1,G,t);
        q1parn = Group_Action_by_Gamma_q(q1par,gam);
        TotalCost(k)=len^2+trapz(linspace(0,1,T),  sum( (q1parn-q2).*(q1parn-q2) )  );
    end
    
    [dmin,indx]=min(TotalCost);
    dsq = dmin;
    slctheta = theta(indx);
    R=AxisAndAngleRotation(u2,slctheta);
    v2R=R*v2;
    w2R=cross(u2,v2R);w2R=w2R/norm(w2R);
    f2=[u2,v2R,w2R];
    A=logm(f2*f1');
    %compute \beta for each \theta
    for j=1:N
        B(:,j)=expm((j-1)*A/(N-1))*u1;
    end
    len=LengthOfTrajectory(B);
    
    %parallel transpot of q2 along \beta to tangent space at u1
    for i=1:T
        q2par(:,i)=BackwardParallelTranslation(B,q2(:,i));
    end
    
    t = linspace(0,1,T);
    [G,T1] = DynamicProgrammingQ2(q1,t,q2par,t,t,t,0);
    gam = interp1(T1,G,t);
    q2parn = Group_Action_by_Gamma_q(q2par,gam);
    epw1 = q2parn - q1;
    
    % % calculate epw1;
    % for i=1:T
    %     q2par(:,i)=BackwardParallelTranslation(B,q2(:,i));
    % end
    % epw1 = q2par - q1;
    dt = 1/T;
    normlized = norm(A*u1,'fro');
    epu1 = len*A*u1/normlized;
    
    % verification
    Rc = zeros(3,1);
    for i=1:T
        Rc = Rc + cross(q1(:,i),epw1(:,i));
    end;
    Rc = Rc*dt;
    Gc = cross(Rc,epu1);
    proj = A*A*u1 - sum(A*A*u1.*u1)*u1;
%     disp(A);

    %plot
    if(disppic == 1)
        figure;clf;
        [x,y,z] = sphere(100);
        h=surf(0.99*x,0.99*y,0.99*z) ;
        axis equal off;
        colormap gray;
        grid off;
        set(h,'LineStyle','none');
        hold on;
        pp1=DensePath(p1,5);
        pp2=DensePath(p2,5);
        plot3(pp1(1,:),pp1(2,:),pp1(3,:),'r','LineWidth',3);
        plot3(pp2(1,:),pp2(2,:),pp2(3,:),'m','LineWidth',3);
        
        for j=1:N
            tt=(j-1)/(N-1);
            qeta=(1-tt)*q1+tt*q2parn;
            for i=1:T
                qgeo(:,i)=ForwardParallelTranslation(B(:,1:j),qeta(:,i));
            end
            pgeo(:,:,j)=q_to_path(qgeo,B(:,j));
        end
        
        geo_0=PiecewiseGeodesic(N,[1 N],[u1,u2]);
        plot3(B(1,:),B(2,:),B(3,:),'y','LineWidth',3);
        plot3(geo_0(1,:),geo_0(2,:),geo_0(3,:),'y--','LineWidth',3);
        
        
        for j=3:2:N-2
            newpp=pgeo(:,:,j);
            plot3(newpp(1,:),newpp(2,:),newpp(3,:),'Color',[1 0 1.0*j/N],'LineWidth',3);
        end
        
    end

end;







