function [p2,q2,path2]=ExponentialMap(p1,epu1,epw1,N)

T=size(p1,2);
%N=21;
disppic = 0;

q1=path_to_q(p1);
u1=p1(:,1); % starting point x
dt = 1/T;

if(norm(epu1)==0)
    p2 = u1;
    q2 = q1 + epw1;
    %trajectory = q_to_path(q2,p2);
    path2 = q_to_path(q2,p2);
else
    %if epu1 != 0
    % epu1: shooting direction for the path x;
    tic
    
    Rc = zeros(3,1);
    for i=1:T
        Rc = Rc + cross(q1(:,i),epw1(:,i));
    end;
    Rc = Rc*dt;
    
    Gc = cross(Rc,epu1);
    
    % This solution is right.
    cnum = u1(1)^2*norm(epu1)^2- Gc(1)*u1(1)-epu1(2)^2-epu1(3)^2;
    elementc = cnum /(epu1(3)*u1(2)-epu1(2)*u1(3));
    elementa = - (epu1(2)-elementc*u1(3))/u1(1);
    elementb = - (epu1(3)+elementc*u1(2))/u1(1);
    
    A = [0 elementa elementb
        -elementa 0 elementc
        -elementb -elementc 0];
    
    %verification;
    for j=1:N
        B(:,j)=expm((j-1)*A/(N-1))*u1;
    end
    
    p2 = expm(A)*u1;
    tmpq2 = q1 + epw1;
    for i=1:T
        q2(:,i)=ForwardParallelTranslation(B,tmpq2(:,i));
    end
    path2 = q_to_path(q2,p2);
    %trajectory = B;
    
    toc
    
    if (disppic == 1)
        
        %plot
        for j=1:N
            tt=(j-1)/(N-1);
            qeta=q1 + tt*epw1;
            for i=1:T
                qgeo(:,i)=ForwardParallelTranslation(B(:,1:j),qeta(:,i));
            end
            pgeo(:,:,j)=q_to_path(qgeo,B(:,j));
        end
        figure(99)
        u2 = expm(A)*u1;
        %plot starting path
        [x,y,z] = sphere(100);clf;
        h=surf(0.96*x,0.96*y,0.96*z) ;
        axis equal off;
        colormap gray;
        grid off;
        set(h,'LineStyle','none');
        hold on;
        plot3(p1(1,:),p1(2,:),p1(3,:),'r','LineWidth',3);
        
        geo_0=PiecewiseGeodesic(N,[1 N],[u1,u2]);
        plot3(B(1,:),B(2,:),B(3,:),'y','LineWidth',3);
        plot3(geo_0(1,:),geo_0(2,:),geo_0(3,:),'y--','LineWidth',3);
        scatter3(B(1,1),B(2,1),B(3,1),120,'or','fill');
        scatter3(B(1,end),B(2,end),B(3,end),120,'sg','fill');
        
        for j=5:5:N
            newpp=pgeo(:,:,j);
            plot3(newpp(1,:),newpp(2,:),newpp(3,:),'Color',[1*(N-j+1)/N 1.0*(j-1)/N 0],'LineWidth',3);
        end
        title('Explicit Result');
    end;
end



