%Example code for simulating and comparing trajectories with phase difference
%author: Zhengwu Zhang
%email: zhengwustat@gmail.com
%date: Jan. 30, 2017

clear all;
close all;

addpath('./RealData')
addpath('./PlotEarth')
addpath('./amplitude_separation')
addpath('./phase_utility_functions')


load simulated_phase_amp.mat;
load gam_sim8;

Ngam = size(gam_sim,1);
figure(10);clf;
plot(gam_sim','linewidt',2);
set(gca,'fontsize',18);


%wrap the data with gamma function, for path 1
path1_s2 = path1{1};
for i=1:Ngam
    path1_wrapped{i} = Group_Action_by_Gamma_p(path1_s2,gam_sim(i,:));
end

%wrap the data with gamma function, for path 2
path2_s2 = path2{1};
for i=1:Ngam
    path2_wrapped{i} = Group_Action_by_Gamma_p(path2_s2,gam_sim(i,:));
end


%resample and smooth the trajectories
for i=1:Ngam
    Sample_N = 50;
    Resample_track = ReSampleSphereTraj(path1_wrapped{i},Sample_N);
    Smoothed_track=SmoothPath(Resample_track,5,0.6);
    Postpreprocess_Path{i} = Smoothed_track;
end;

for i=1:Ngam
    Sample_N = 50;
    Resample_track = ReSampleSphereTraj(path2_wrapped{i},Sample_N);
    Smoothed_track=SmoothPath(Resample_track,5,0.6);
    Postpreprocess_Path{Ngam+i} = Smoothed_track;
end;


%plot all paths on sphere
figure(1);clf;
[x,y,z] = sphere(100);
h=surf(x,y,z) ;
axis equal off;
colormap gray;
grid off;
set(h,'LineStyle','none');

figure(2);clf;
for j=1:2*Ngam
    figure(1);
    hold on;
    p1 = Postpreprocess_Path{j};
    pp1=DensePath(p1,5);
    color = 'r';
    if(j>Ngam)
        color = 'b';
    end
    plot3(pp1(1,:),pp1(2,:),pp1(3,:),color,'LineWidth',3);
    
    figure(2);
    hold on;
    f1 = path_to_function(p1);
    plot(1:T,f1,color,'LineWidth',2);
end

%%%%%%%%%%%%%%%%%%% simulate the phase %%%%%%%%%%%%%%%%%%%%
% simulate the phase part with random gam function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%consider the phase difference
c1 = 0;  %weight for amp component
c2 = 1; % weight for pha component
dT = 1/(1-T);

Dmat_all = zeros(2*Ngam,2*Ngam);
gamID = linspace(0,1,T);
psiid = sqrt(diff(gamID)/dT);


figure(1);clf;
for i=1:2*Ngam
    path1 = Postpreprocess_Path{i};
    for j=i+1:2*Ngam
        path2 = Postpreprocess_Path{j};
        [dmin_amp,~,gamI]=GeodesicsWithRegestration_coordecent(path1,path2,1);
        psi = sqrt(diff(gamI)/dT);
        d_phase =  acos(sum(psiid.*psi)*dT);
        
        Dmat_amp(i,j) = dmin_amp;
        Dmat_phase(i,j) = d_phase;
        
        Dmat_all(i,j) = c1*dmin_amp + c2*d_phase;
        figure(1),hold on;
        plot(linspace(0,1,T),gamI);
        %title(sprintf('(%d,%d) d = %0.2f',i,j,dmin));
        %pause(0.5);
    end
end

Dmat_amp(2*Ngam,2*Ngam) = 0;
Dmat_phase(2*Ngam,2*Ngam) = 0;

c1 = 0.5;  %weight for amp.
c2 = 0.5; % weight for phase

Dmat_all = c1*Dmat_amp + c2*Dmat_phase;
figure, imagesc(Dmat_all+Dmat_all');
set(gca,'fontsize',18);
colorbar;

