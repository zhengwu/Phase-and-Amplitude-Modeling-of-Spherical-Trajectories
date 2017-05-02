function [p2,q2,trajectory]=ExponentialMap_Numerical(p1,epu1,epw1,N)
% function for numerical exponential map
% input: p1 -- starting point p1 (a path)
%      : (epu1, epw1) shooting vectors. 
T=size(p1,2);
%N=21;
disppic = 1;

q1=path_to_q(p1);
u1=p1(:,1); % starting point x
dt = 1/T;

if(norm(epu1)==0)
    p2 = u1;
    q2 = q1 + epw1;
    trajectory = q_to_path(q2,p2);
else
   % numerical constructing the exponential map. 
   tic
    epsilon = 1/N;
    B(:,1) = u1;
    B(:,2) = Exp_Sphere(B(:,1),epsilon*epu1); % exponential map;
    B(:,2) = B(:,2)/norm(B(:,2));
    qgeo(:,:,1) = q1;
    for j=1:T
        vparall(:,j) = ForwardParallelTranslation(B(:,1:2),q1(:,j));
        wparall(:,j) = ForwardParallelTranslation(B(:,1:2),epw1(:,j));
    end;
    qgeo(:,:,2) = vparall + epsilon*wparall;
    xs(:,1) = epu1;
    
    Rc = zeros(3,1);
    for j=1:T
        Rc = Rc + cross(q1(:,j),epw1(:,j));
    end;
    Rc = Rc*dt;
    Gc = cross(Rc,xs(:,1));
    dxs(:,1) = Gc;
    
    for i=2:N-1
        xs(:,i) = ForwardParallelTranslation(B(:,i-1:i),xs(:,i-1)+epsilon*dxs(:,i-1));
        %B(:,i+1) = Exp_Sphere(B(:,i),epsilon*xs(:,i-1)); % baseline
        B(:,i+1) = Exp_Sphere(B(:,i),epsilon*xs(:,i)); % baseline
        B(:,i+1) = B(:,i+1)/norm(B(:,i+1));
        for j=1:T % update geodesic; 
            vparall(:,j) = ForwardParallelTranslation(B(:,1:i+1),q1(:,j));
            wparall(:,j) = ForwardParallelTranslation(B(:,1:i+1),epw1(:,j));
        end;
        qgeo(:,:,i+1) =  vparall + epsilon*(i+1)*wparall;% path on tangent space.
        
        % calculate the dxs;
        for j=1:T
            vparall(:,j) = ForwardParallelTranslation(B(:,1:i),q1(:,j));
            wparall(:,j) = ForwardParallelTranslation(B(:,1:i),epw1(:,j));
        end;
        Rc = zeros(3,1);
        for j=1:T
            Rc = Rc + cross(vparall(:,j),wparall(:,j));
        end;
        Rc = Rc*dt;
        Gc = cross(Rc,xs(:,i));
        dxs(:,i) = Gc;
    end;
  end;
  
   q2 = qgeo(:,:,end);
   p2 = q_to_path(qgeo(:,:,end),B(:,end));
   trajectory = B;
   
   toc
    if (disppic == 1)
        
        %plot
        for j=1:N
            pgeo(:,:,j)=q_to_path(qgeo(:,:,j),B(:,j));
        end
        figure(100)
        %plot starting path
        [x,y,z] = sphere(100);clf;
        h=surf(0.96*x,0.96*y,0.96*z) ;
        axis equal off;
        colormap gray;
        grid off;
        set(h,'LineStyle','none');
        hold on;
        pp1=DensePath(p1,5);
        plot3(pp1(1,:),pp1(2,:),pp1(3,:),'r','LineWidth',3);
        
        geo_0=PiecewiseGeodesic(N,[1 N],[u1,B(:,end)]);
        plot3(B(1,:),B(2,:),B(3,:),'y','LineWidth',3);
        plot3(geo_0(1,:),geo_0(2,:),geo_0(3,:),'y--','LineWidth',3);
        
        scatter3(B(1,1),B(2,1),B(3,1),120,'or','fill');
        scatter3(B(1,end),B(2,end),B(3,end),120,'sg','fill');
        
        for j=5:5:N
            newpp=pgeo(:,:,j);
            plot3(newpp(1,:),newpp(2,:),newpp(3,:),'Color',[1*(N-j+1)/N 1.0*(j-1)/N 0],'LineWidth',3);
        end
        title('Numerical Result')
    end;
end



