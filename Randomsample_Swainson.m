%Example code for fitting a mixture model and sampling
% from the fitted distribution
%author: Zhengwu Zhang
%email: zhengwustat@gmail.com
%date: Jan. 30, 2017

clear all;
close all;

addpath('./RealData')
addpath('./PlotEarth')
addpath('./amplitude_separation')
addpath('./phase_utility_functions')

load swainson.mat;
TotalN = length(swainson_date);

%Plot the earth and trajectories
figure(1);hold on;
globe([],'earth_1600.png');
arg='100,''yellow'',''fill'',''markeredgecolor'',''black'' ';
hold on;

T = 100;
N = 10;
for i=1:N
    Resampled_track = ReSampleSphereTraj( swainson_cord{3*i+3},T);
    Smoothed_track=SmoothPath(Resampled_track,7,1);
    path{i} = Smoothed_track;
end

for i=1:N
    hold on;
    color = 'yellow';
    tmppath = path{i};
    plot3(tmppath(1,:),tmppath(2,:),tmppath(3,:),color,'linewidth',0.6);
end;

%%%%%%random sample based on the seperate modeling 
%(of the two parts of shooting vectors) strategy
%[pcapath,mup,muq,mupath]= amples_covm_pca_randsamples_separately(path,'slow');



%%%%%%random sampling based on the joint modeling 
%(of the two parts of shooting vectors) strategy
[pcapath,mup,muq,mupath]= amples_covm_pca_randsamples_jointly(path,'slow');



