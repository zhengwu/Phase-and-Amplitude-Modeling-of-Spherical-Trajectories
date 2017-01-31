%Example code for calculating the Karcher mean of given a set of
%trajectories
%author: Zhengwu Zhang
%email: zhengwustat@gmail.com
%date: Jan. 30, 2017

clear all;
close all;

addpath('./RealData')
addpath('./PlotEarth')
addpath('./amplitude_separation')

load swainson.mat;
TotalN = length(swainson_date);

%Plot Earth
figure(100);hold on;
globe([],'earth_1600.png');
arg='100,''yellow'',''fill'',''markeredgecolor'',''black'' ';
hold on;

%Resample the data
T = 100;
N = 10;
for i=1:N
    Resampled_track = ReSampleSphereTraj( swainson_cord{3*i+3},T);
    Smoothed_track=SmoothPath(Resampled_track,7,1);
    path{i} = Smoothed_track;
end

%Plot the trajectories
for i=1:N
    hold on;
    color = 'yellow';
    tmppath = path{i};
    plot3(tmppath(1,:),tmppath(2,:),tmppath(3,:),color,'linewidth',0.6);
end;

%pause to show figure 100;
pause(1);

%Calculate the Euclidean mean
average_path = Euclidean_Mean(path, 1);

%Calculate the Karcher mean
[mup,muq,mupath]= KarcherMean(path,'slow');

figure(101);hold on;
globe([],'earth_1600.png');
arg='100,''yellow'',''fill'',''markeredgecolor'',''black'' ';
hold on;
plot3(mupath(1,:),mupath(2,:),mupath(3,:),'m','LineWidth',3);

for i=1:N
     [pathn{i},indx,gam]=Allignp1top2(path{i},mupath);
end;

%%%%Figure 102 plots the Karcher mean and 
% the cross-sectional variance at each point along the Karcher mean
figure(102);clf;hold on;
globe([],'earth_1600.png');
arg='100,''yellow'',''fill'',''markeredgecolor'',''black'' ';
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
    Total_var_align(i) = trace(K);
    if mod(i,10)==0 || i==1
        lam1=sqrt(S(1,1));
        lam2=sqrt(S(2,2));
        the=linspace(0,2*pi,100);
        [xthe]=1*lam1*sin(the);
        [ythe]=1*lam2*cos(the);
        yyy=U*[xthe;ythe;zeros(1,100)]+repmat(mupath(:,i),1,100);
        plot3(yyy(1,:),yyy(2,:),yyy(3,:),'y','LineWidth',3);
    end
end;

