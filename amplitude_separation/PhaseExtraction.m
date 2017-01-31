function [GAM]= PhaseExtraction(mupath,patharray)
% function for calculating Karcher Mean of given trajectories on S2
% Input -- patharray: cell type, each cell contain a trajectory with dim = d*T, 
%                 d is the dimension of the path (usually d = 3),
%                 T is the # of sampling points on one trajectory

%Output -- mup: the starting point of the mean trajectory
%       -- muq: the TSRVF of the mean trajectory
%       -- mupath: The path of the mean trajectory
N = length(patharray);
[d,T] = size(patharray{1});

%transform path to TSRVF
for i=1:N
    cpath = patharray{i};
    ptarray{i} = cpath(:,1);
    qarray{i} = path_to_q(cpath);
end;

mup = mupath(:,1);
muq = path_to_q(mupath);

for i=1:N
    [epmu(:,i),epmw(:,:,i),dsq,gam]=InverseExponentialMapPhase(patharray{i},mupath);
     GAM(i,:) = gam;
end;
