%Example code for calculating the Karcher mean of given a set of
%trajectories
%author: Zhengwu Zhang
%email: zhengwustat@gmail.com
%date: Jan. 30, 2017


clear all;
close all;

load Hurricane_subset1.mat;
TotalN = length(Hurri_Data);

% resample the trajectories
T = 100;
N = 7;
for i=1:N
    Resampled_track = ReSampleSphereTraj(Hurri_Data{i},T);
    Smoothed_track=SmoothPath(Resampled_track,7,1);
    path{i} = Smoothed_track;
end

%calculate the Euclidean mean
average_path = Euclidean_Mean(path, 1);
%calculate the Karcher Mean
[mup,muq,mupath]= KarcherMean(path,'slow');

figure(100);hold on;
globe([],'earth_1600.png');
arg='100,''yellow'',''fill'',''markeredgecolor'',''black'' ';
hold on;
plot3(mupath(1,:),mupath(2,:),mupath(3,:),'g','LineWidth',3);
for i=1:N
    hold on;
    color = 'white';
    tmppath = path{i};
    plot3(tmppath(1,:),tmppath(2,:),tmppath(3,:),color,'linewidth',3);
end;


for i=1:N
     [pathn{i},indx,gam]=Allignp1top2(path{i},mupath);
end;
figure(111);clf;hold on;
globe([],'earth_1600.png');
arg='100,''yellow'',''fill'',''markeredgecolor'',''black'' ';
hold on;
hold on;
PLOT(mupath,1,'w',3);

% calculate the convariance and plot the convariance
for i=1:T
    muX=mupath(:,i);
    for j=1:N
        Dat(j,:)=InverseExp_Sphere(muX,pathn{j}(:,i));
    end
    K=cov(Dat);
    [U,S,V]=svd(K);
    Total_var_align(i) = trace(S);
    if mod(i,15)==0 || i==1
        lam1=sqrt(S(1,1));
        lam2=sqrt(S(2,2));
        the=linspace(0,2*pi,100);
        [xthe]=1*lam1*sin(the);
        [ythe]=1*lam2*cos(the);
        yyy=U*[xthe;ythe;zeros(1,100)]+repmat(mupath(:,i),1,100);
        plot3(yyy(1,:),yyy(2,:),yyy(3,:),'y','LineWidth',3);
    end
end;

