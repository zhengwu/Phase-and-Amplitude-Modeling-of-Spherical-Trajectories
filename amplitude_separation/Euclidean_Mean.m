function average_path = Euclidean_Mean(path, disp)

N = length(path);
[d,T] = size(path{1});

%average
average_path = zeros(3,T);
for i=1:N
    average_path = average_path + path{i};
end;
average_path = average_path/N;
for i=1:T
    average_path(:,i) = average_path(:,i)/norm(average_path(:,i));
end;

if(disp == 1)
    %plot the data and convariance matrix
    figure(101);hold on;
    globe([],'earth_1600.png');
    arg='100,''yellow'',''fill'',''markeredgecolor'',''black'' ';
    hold on;
    hold on;
    PLOT(average_path,1,'w',3);
    
    % calculate the convariance and plot the convariance
    for i=1:T
        muX=average_path(:,i);
        for j=1:N
            Dat(j,:)=InverseExp_Sphere(muX,path{j}(:,i));
        end
        K=cov(Dat);
        [U,S,V]=svd(K);
        
        Total_Var(i) = trace(S);
        if mod(i,15)==0 || i==1
            lam1=sqrt(S(1,1));
            lam2=sqrt(S(2,2));
            the=linspace(0,2*pi,100);
            [xthe]=1*lam1*sin(the);
            [ythe]=1*lam2*cos(the);
            yyy=U*[xthe;ythe;zeros(1,100)]+repmat(average_path(:,i),1,100);
            plot3(yyy(1,:),yyy(2,:),yyy(3,:),'y','LineWidth',3);
        end
    end;
end;
figure(102);clf;
plot(1:T,Total_Var,'--r','linewidth',2.5);
