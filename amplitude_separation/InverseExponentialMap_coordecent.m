function [epu1,epw1,dsq]=InverseExponentialMap_coordecent(p1,p2)

T=size(p1,2);
N=21;
disppic = 0;

q1=path_to_q(p1);
q2=path_to_q(p2);


u1=p1(:,1);
u2=p2(:,1);


if ( norm(u1-u2)<0.001 )
    
    % The starting points are the same
    %compute the cost
    len=0;
    t = linspace(0,1,T);
    [G,T1] = DynamicProgrammingQ2(q2,t,q1,t,t,t,0);
    gam = interp1(T1,G,t);
    gamI = invertGamma(gam);
    gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
    p2n = Group_Action_by_Gamma_p(p2,gamI);
    q2n = path_to_q(p2n);
    dsq=len^2+trapz(linspace(0,1,T),  sum( (q1-q2n).*(q1-q2n)));
    epw1 = q2n - q1;
    epu1 = zeros(3,1);
else 
    
    v1=InverseExp_Sphere(u1,p1(:,8));
    v1=v1/norm(v1);
    v2=InverseExp_Sphere(u2,p2(:,8));
    v2=v2/norm(v2);
    
    w1=cross(u1,v1);
    w1=w1/norm(w1);
    w2=cross(u2,v2);
    w2=w2/norm(w2);
    
    f1=[u1,v1,w1];
    f2=[u2,v2,w2];
    
    %inital theta
    M=5;
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
        gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
        p2n = Group_Action_by_Gamma_p(p2,gamI);
        q2n = path_to_q(p2n);
        TotalCost(k)=len^2+trapz(linspace(0,1,T),  sum( (q1par-q2n).*(q1par-q2n) )  );
    end
    [dmin,indx]=min(TotalCost);
    
    % Steepest-Descent for searching the optimal theta
    h = 0.04;
    delta = 1;
    epsilon = 0.002;
    xnew = theta(indx); %inital
    maxit = 60;
    iter = 0;fp = 10;
    while(abs(fp)>epsilon)
        iter = iter + 1;
        if(iter>maxit)
            break;
            display('Max interation has reached');
        end;
        
        x = xnew;
        xh = x+h;
        
        R = AxisAndAngleRotation(u2,x);
        Rh = AxisAndAngleRotation(u2,xh);
        
        v2R = R*v2;
        v2Rh = Rh*v2;
        
        w2R=cross(u2,v2R);w2R=w2R/norm(w2R);
        w2Rh=cross(u2,v2Rh);w2Rh=w2Rh/norm(w2Rh);
        
        f2=[u2,v2R,w2R];
        f2h = [u2,v2Rh,w2Rh];
        
        A=logm(f2*f1');
        Ah = logm(f2h*f1');
        
        %compute \beta for both
        for j=1:N
            B(:,j)=expm((j-1)*A/(N-1))*u1;
        end
        for j=1:N
            Bh(:,j)=expm((j-1)*Ah/(N-1))*u1;
        end
        
        for i=1:T
            q1par(:,i)=ForwardParallelTranslation(B,q1(:,i));
        end
        
        for i=1:T
            q1parh(:,i)=ForwardParallelTranslation(Bh,q1(:,i));
        end
        
        %compute the cost
        len=LengthOfTrajectory(B);
        lenh=LengthOfTrajectory(Bh);
        t = linspace(0,1,T);
        
        [G,T1] = DynamicProgrammingQ2(q2,t,q1par,t,t,t,0);
        gam = interp1(T1,G,t);
        gamI = invertGamma(gam);
        gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
        p2n = Group_Action_by_Gamma_p(p2,gamI);
        q2n = path_to_q(p2n);
        
         
        [G,T1] = DynamicProgrammingQ2(q2,t,q1parh,t,t,t,0);
        gam = interp1(T1,G,t);
        gamI = invertGamma(gam);
        gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
        p2nh = Group_Action_by_Gamma_p(p2,gamI);
        q2nh = path_to_q(p2nh);
        
        TotalCost=len^2+trapz(linspace(0,1,T),  sum( (q1par-q2n).^2 )  );
        TotalCosth=lenh^2+trapz(linspace(0,1,T),  sum( (q1parh-q2nh).^2 )  );
        
        fp = (TotalCosth - TotalCost)/h;
        
        xnew = x - delta*fp;
        xtheta(iter) = x;
        Cvalue(iter) = TotalCost;
        dvalue(iter) = fp;
    end
    minCost = TotalCosth;
    slctheta = xnew;
    dsq = TotalCosth;
    
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
    
    %parallel transpot of q1 along \beta to tangent space at u2
    for i=1:T
        q2par(:,i)=BackwardParallelTranslation(B,q2(:,i));
    end
    t = linspace(0,1,T);
    [G,T1] = DynamicProgrammingQ2(q1,t,q2par,t,t,t,0);
    gam = interp1(T1,G,t);
    gamI = invertGamma(gam);
    gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
    p1n = Group_Action_by_Gamma_p(p1,gamI);
    q1n = path_to_q(p1n); 
   
    epw1 = q2par - q1n;
    
    dt = 1/T;
    normlized = norm(A*u1,'fro');
    epu1 = len*A*u1/normlized;
    
    % verification
    Rc = zeros(3,1);
    for i=1:T
        Rc = Rc + cross(q1(:,i),epw1(:,i));
    end;
    Rc = Rc*dt;
    Gc = cross(Rc,epu1)
    proj = A*A*u1 - sum(A*A*u1.*u1)*u1
    disp(A);
    
    
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
            qeta=(1-tt)*q1+tt*q2npar;
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




