function [mup,muq,mupath]= KarcherMean(patharray,str)
% function for calculating Karcher Mean of given trajectories on S2
% Input -- patharray: cell type, each cell contain a trajectory with dim = d*T, 
%                 d is the dimension of the path (usually d = 3),
%                 T is the # of sampling points on one trajectory

%Output -- mup: the starting point of the mean trajectory
%       -- muq: the TSRVF of the mean trajectory
%       -- mupath: The path of the mean trajectory
N = length(patharray);
[d,T] = size(patharray{1});

%parameters for the whole algorithm
maxiter = 25; % max iterations for calculating the Karcher Mean;
thrdv = 0.02;
thrdd = 0.02;
eps = 0.5;
plotfig=0;
iter = 2;

%transform path to TSRVF
for i=1:N
    cpath = patharray{i};
    ptarray{i} = cpath(:,1);
    qarray{i} = path_to_q(cpath);
end;

% Initial mu as one of the paths;
mupath = patharray{1};
mup = mupath(:,1);
muq = path_to_q(mupath);

while (iter<maxiter)
    iter = iter + 1
    % calculate the inverse expenential map;
    epmu = zeros(d,N);
    epmw = zeros(d,T);
    tmpsumd = 0;
    for i=1:N
        if(str =='fast')
            [epmu(:,i),epmw(:,:,i),dsq]=InverseExponentialMapfast(mupath,patharray{i});
        else
            [epmu(:,i),epmw(:,:,i),dsq]=InverseExponentialMap(mupath,patharray{i});        
        end;
        disp(dsq);
        tmpsumd = tmpsumd + dsq;
    end;

    epm_ubar = sum(epmu,2)/N;
    epm_wbar = sum(epmw,3)/N;
    
    sumd(iter) = tmpsumd;
    norm_ubar(iter) = norm(epm_ubar);
    norm_wbar(iter) = trapz(linspace(0,1,T),  sum( epm_wbar.^2 )  );
    
    norm_sumv = norm_ubar(iter) + norm_wbar(iter);
    summd_cum = abs(sumd(iter)-sumd(iter-1)) %+ abs(sumd(iter-1)-sumd(iter-2));
    if (norm_sumv>thrdv) && (summd_cum>thrdd)
        %update mu
         update_u = eps*epm_ubar;
         update_w = eps*epm_wbar;
         [mup,muq,mupath] = ExponentialMap(mupath,update_u,update_w,21);
         
         %figure(100);hold on;
         %plot3(mupath(1,:),mupath(2,:),mupath(3,:),'m','LineWidth',3);
    else
        break
    end
    
end;

if(plotfig ==1)
figure, plot(2:iter,norm_ubar(2:iter),'g','linewidth',2);
title('Plot of norm ubar')
figure, plot(2:iter,norm_wbar(2:iter),'linewidth',2);
title('Plot of norm wbar')
figure, plot(2:iter,sumd(2:iter),'r','linewidth',2);
title('Plot of norm sumd2')
end

