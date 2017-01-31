%Example code for calculating the distance matrix for given hurricanes trajectories 
%author: Zhengwu Zhang
%email: zhengwustat@gmail.com
%date: Jan. 30, 2017

clear all;
close all;

addpath('./RealData')
addpath('./PlotEarth')
addpath('./amplitude_separation')

%myPool = parpool();

%load data
% load hurricane_after1969_coordinate;
% TotalN = length(hurr_af1969_year);

load hurricane_20n_35n_coordinate;
TotalN = length(hurr_20n_35n_year);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----Translate The Coordinate From Spherical(theta,phi) to Cartesian --%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=1.02;
for i=1:TotalN
%     tmp_theta = hurr_af1969_theta{i};
%     tmp_phi = hurr_af1969_phi{i};

    tmp_theta = hurr_20n_35n_theta{i};
    tmp_phi = hurr_20n_35n_phi{i};

    % sample points along the trajectory
    N = length(tmp_theta);
   
    [X(1,:),X(2,:),X(3,:)] =s2c(tmp_theta,tmp_phi);
    X=r*X;
    hurri_path_cartesian{i} = X;
    clear X;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- Pre-process The Data--%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resample and smooth the track;
for i=1:TotalN
    Sample_N = 50;
    Resample_track = ReSampleSphereTraj(hurri_path_cartesian{i},Sample_N);
    Smoothed_track=SmoothPath(Resample_track,5,0.6);
    Postpreprocess_Hurri_Path{i} = Smoothed_track;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- Calculate the Distance Matrix Using OUR Method--%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DistM = zeros(TotalN,TotalN);
for i=1:TotalN
    parfor j=i+1:TotalN
        track1 = Postpreprocess_Hurri_Path{i};
        track2 = Postpreprocess_Hurri_Path{j};
%         tic
%         [dminfastg,ind1,gamI]=GeodesicsWithRegestration_gridsearch(track1,track2,1);
%         toc
        tic
        [dminfastc,ind1,gamI]=GeodesicsWithRegestration_coordecent(track1,track2,1);
        toc
        DistM(i,j) = dminfastc;
    end;
    
end;
DistM = DistM + DistM';

save DistM_Hurricane DistM;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- Calculate the distance matrix using the method in Su et al. 2014 --%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = [1,0,0]';
T = Sample_N;
DistMSU = zeros(TotalN,TotalN);
dim=3;
for i=1:TotalN
    parfor j=i+1:TotalN
        track1 = Postpreprocess_Hurri_Path{i};
        track2 = Postpreprocess_Hurri_Path{j};
        
         q1=path_to_q_su(track1,c);
         q2=path_to_q_su(track2,c);
         
         t = linspace(0,1,T);
         [G,T1] = DynamicProgrammingQ2(q2,t,q1,t,t,t,0);
         gam = interp1(T1,G,t);
         gamI = invertGamma(gam);
         gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
         gamI_dev = gradient(gamI, 1/(T-1));
         q2n = interp1(t', q2', (t(end)-t(1)).*gamI' + t(1), 'linear', 'extrap')'.*(ones(dim,1)*sqrt(gamI_dev));
         
         dmin=sqrt( sum(sum( (q1-q2n).*(q1-q2n))/(T-1) ));
         DistMSU(i,j) = dmin;
        
    end;
    
end;
DistMSU = DistMSU + DistMSU';

save DistMSu_Hurricane DistMSU;
