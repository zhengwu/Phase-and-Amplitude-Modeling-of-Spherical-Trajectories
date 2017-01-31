%Example code for comparing two trajectories
%author: Zhengwu Zhang
%email: zhengwustat@gmail.com
%date: Jan. 30, 2017

clear all;
close all;


addpath('./RealData')
addpath('./PlotEarth')
addpath('./amplitude_separation')

%load data;
load swainson.mat;
TotalN = length(swainson_date);

%two tracks
track1 =  swainson_cord{4};
track2 =  swainson_cord{10};

% comparing these two tracks;
n1 = size(track1,2);
n2 = size(track2,2);

% Resample the track;
N = 50;
X1 = ReSampleSphereTraj(track1,N);
X2 = ReSampleSphereTraj(track2,N);

% smooth the track;
smoothed_X1=SmoothPath(X1,5,0.6);
smoothed_X2=SmoothPath(X2,5,0.6);


globe([],'BlueMarble.png');
arg='100,''yellow'',''fill'',''markeredgecolor'',''black'' ';
hold on;
hold on,plot3(smoothed_X1(1,:),smoothed_X1(2,:),smoothed_X1(3,:),'g','LineWidth',3);
hold on,plot3(smoothed_X2(1,:),smoothed_X2(2,:),smoothed_X2(3,:),'g','LineWidth',3);


% compare two trajectories without registration
[dmin5,ind5]=GeodesicsWithoutRegestration(X1,X2,6);

% compare two trajectories with registration
[dminfast,ind1,gamI]=GeodesicsWithRegestration_gridsearch(smoothed_X1,smoothed_X2,1);
%[dminfast,ind1,gamI]=GeodesicsWithRegestration_coordecent(X1,X2,1);

