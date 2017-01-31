%Example code for phase and amplitude seperation 
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


T = 100;
N = 10;
for i=1:N
    Resampled_track = ReSampleSphereTraj( swainson_cord{3*i+3},T);
    Smoothed_track=SmoothPath(Resampled_track,7,1);
    path{i} = Smoothed_track;
end

%Plot the earth and trajectories
figure(1);hold on;
globe([],'earth_1600.png');
arg='100,''yellow'',''fill'',''markeredgecolor'',''black'' ';
hold on;
for i=1:N
    hold on;
    color = 'yellow';
    tmppath = path{i};
    plot3(tmppath(1,:),tmppath(2,:),tmppath(3,:),color,'linewidth',0.6);
end;

%calculate the Karcher mean
[mup,muq,mupath]= KarcherMean(path,'slow');

%amplitude component - Amp_comp
for i=1:N
     [Amp_comp{i},indx,gam]=Allignp1top2(path{i},mupath);
end;

[Phase_comp]= PhaseExtraction(mupath,path);
figure(2);clf;
T=100;
t = (0:T-1)/(T-1);
plot(t,Phase_comp','linewidth',2);
set(gca,'fontsize',22);