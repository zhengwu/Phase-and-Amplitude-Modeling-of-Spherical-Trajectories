function [final_dist,indx,gamI]=GeodesicsWithRegestration_gridsearch(p1,p2,fig_id)

T=size(p1,2);
N=21; % sample points of baseline
M=60; % grid search parameter; 60 search point between (0,2*pi);

plotfig = 1;

if (plotfig>1)
    % plot the original path, p1 and p2
    figure(fig_id);clf;
    [x,y,z] = sphere(100);
    h=surf(0.99*x,0.99*y,0.99*z) ;
    axis equal off;
    colormap gray;
    grid off;
    set(h,'LineStyle','none');
    hold on;
    plot3(p1(1,:),p1(2,:),p1(3,:),'r','LineWidth',3);
    plot3(p2(1,:),p2(2,:),p2(3,:),'m','LineWidth',3);
end;

q1=path_to_q(p1);
q2=path_to_q(p2);


u1=p1(:,1);
u2=p2(:,1);

v1=InverseExp_Sphere(u1,p1(:,10));
v1=v1/norm(v1);
v2=InverseExp_Sphere(u2,p2(:,10));
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
    gamI = invertGamma(gam);
    gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
    p2nh = Group_Action_by_Gamma_p(p2,gamI);
    q2nh = path_to_q(p2nh);

    %     q2n = Group_Action_by_Gamma(q2,gamI);
    final_dist=sqrt(trapz(linspace(0,1,T),  sum( (q1-q2nh).*(q1-q2nh) ) ));
    indx = 0;
else  
   
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
        % apply gamI to path 2
        gamI = invertGamma(gam);
        gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
        p2nh = Group_Action_by_Gamma_p(p2,gamI);
        q2nh = path_to_q(p2nh);
        
        TotalCost(k)=len^2+trapz(linspace(0,1,T),  sum( (q1par-q2nh).*(q1par-q2nh) )  );
        Part1(k) = len^2;
        Part2(k) = trapz(linspace(0,1,T),  sum( (q1par-q2nh).*(q1par-q2nh) )  );
    end
    
    
if (plotfig>1)
    % plot the grid search result: cost function w.r.p to the search
    % parameter
    figure(110);clf;
    plot(TotalCost);
    hold on,plot(Part1,'-r','linewidth',3)
    hold on, plot(Part2,'-g','linewidth',3)
    title('After Registration')
end;
    
    % Find the optimal \beta
    [dmin,indx]=min(TotalCost);
    R=AxisAndAngleRotation(u2,theta(indx));
    v2R=R*v2;
    w2R=cross(u2,v2R);w2R=w2R/norm(w2R);
    f2=[u2,v2R,w2R];
    A=logm(f2*f1');
    
    %compute the baseline
    for j=1:N
        B(:,j)=expm((j-1)*A/(N-1))*u1;
    end
    
    %parallel transpot of q1 along \beta to tangent space at u2
    for i=1:T
        q1par(:,i)=ForwardParallelTranslation(B,q1(:,i));
    end
    
    [G,T1] = DynamicProgrammingQ2(q2,t,q1par,t,t,t,0);
    gam = interp1(T1,G,t);
    gamI = invertGamma(gam);
    gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
    p2n = Group_Action_by_Gamma_p(p2,gamI);
    q2n = path_to_q(p2n);
    
    len=LengthOfTrajectory(B);
    final_dist = sqrt(len^2 + trapz(linspace(0,1,T),  sum( (q1par-q2n).^2 )));
    
    for j=1:N
        tt=(j-1)/(N-1);
        qeta=(1-tt)*q1par+tt*q2nh;
        for i=1:T
            qgeo(:,i)=ForwardParallelTranslation(B(:,end:-1:j),qeta(:,i));
        end
        pgeo(:,:,j)=q_to_path(qgeo,B(:,j));
    end
    
if (plotfig>0)
    figure(fig_id);
    [x,y,z] = sphere(100);
    h=surf(0.99*x,0.99*y,0.99*z) ;
    axis equal off;
    colormap gray;
    grid off;
    set(h,'LineStyle','none');
    hold on;
    geo_0=PiecewiseGeodesic(N,[1 N],[u1,u2]);
    plot3(B(1,:),B(2,:),B(3,:),'y','LineWidth',3);
    hold on;
    plot3(geo_0(1,:),geo_0(2,:),geo_0(3,:),'y--','LineWidth',3);
    for j=6:5:N
        newpp=pgeo(:,:,j);
        plot3(newpp(1,:),newpp(2,:),newpp(3,:),'Color',[1 0 1.0*j/N],'LineWidth',3);
        hold on;
    end
    figure(fig_id+1)
    plot(gamI,'linewidth',2);
end
    
end;





