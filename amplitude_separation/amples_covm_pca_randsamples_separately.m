function [pcapath,mup,muq,mupath]= amples_covm_pca_randsamples_separately(patharray,str)
% function for calculating covariance matrices of given trajectories on S2
% we treat two components seperately in this function
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

iter = 2;

%transform path to TSRVF
for i=1:N
    cpath = patharray{i};
    ptarray{i} = cpath(:,1);
    qarray{i} = path_to_q(cpath);
end;

% Initialize mu as one of the paths;
mupath = patharray{1};
mup = mupath(:,1);
muq = path_to_q(mupath);

while (iter<maxiter)
    iter = iter + 1
    % calculate the inverse exponential map;
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
         
         figure(100);hold on;
         plot3(mupath(1,:),mupath(2,:),mupath(3,:),'m','LineWidth',3);
    else
        break
    end
    
end;

keyboard
%calculate the shooting vectors
for i=1:N
    if(str =='fast')
        [epmu(:,i),epmw(:,:,i),dsq]=InverseExponentialMapfast(mupath,patharray{i});
    else
        [epmu(:,i),epmw(:,:,i),dsq]=InverseExponentialMap(mupath,patharray{i});
    end;
end;

%form a new coordinate system on T_mu_p(S^2)
v1 = epmu(:,2);
v1 = v1/norm(v1);
if(sum(mup.*v1)>0.01)
    disp('Tangent vectors do not perpendicular to vector mu...')
    keyboard
end;

v2 = epmu(:,6)+rand(3,1);
v2 = v2 - mup*dot(v2,mup); %make v2 prependicular to v1;
v2 = v2/norm(v2);
v2 = v2 - v1*dot(v1,v2);
v2 = v2/norm(v2);

%represent shooting vectors using new coordinat system
for i=1:N
    new_epmu(1,i) = dot(epmu(:,i),v1);
    new_epmu(2,i) = dot(epmu(:,i),v2);
    
    for j=1:T
        new_epmw(1,j,i) = dot(epmw(:,j,i),v1);
        new_epmw(2,j,i) = dot(epmw(:,j,i),v2);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% perform PCA analysis for epmu and epmw jointly 
%PCA analysis for each shooting components jointly 
%first component;
Sp = zeros(2,2);
mu_new_epmu = mean(new_epmu,2);
for i=1:N
    tmp = new_epmu(:,i)-mu_new_epmu;
    Sp = Sp + tmp*tmp';
end;
Sp = Sp/(N-1);

[Up,Sigmau]=svd(Sp);

%second component;
for i=1:N
vnew_epmw(:,i) = reshape(new_epmw(:,:,i)',[2*T,1]); %firt 1:T is x, second 1:T is y;
end
mu_vnew_epmw = zeros(size(mean(vnew_epmw,2)));%mean(vnew_epmw,2);
mu_vnew_epmw = mean(vnew_epmw,2);
K = cov(vnew_epmw');
[Um,Sigmam]=svd(K);

%display the PCA
tau = -1:(1/3):1;
pcu = 1;
pcw = 1;
ind_u = 1;
ind_w = 0;

for i=1:length(tau)
      new_pca_epmu(:,i) = ind_u*(mu_new_epmu + tau(i)*sqrt(Sigmau(pcu,pcu))*Up(:,pcu));
      new_pca_epmw(:,:,i) = ind_w*(reshape(mu_vnew_epmw + tau(i)*sqrt(Sigmam(pcw,pcw))*Um(:,pcw),[T,2])');
end

%represent the PCA using old coordinate system
for i=1:length(tau)
    pca_epmu(:,i) = v1*new_pca_epmu(1,i) + v2*new_pca_epmu(2,i);
    
    for j=1:T
        pca_epmw(:,j,i) = v1*new_pca_epmw(1,j,i) + v2*new_pca_epmw(2,j,i);
    end;
end

%reconstruct curves
for i=1:length(tau)
    [pcap(:,i),pcaq(:,:,i),pcapath(:,:,i)] = ExponentialMap(mupath,pca_epmu(:,i),pca_epmw(:,:,i),21);
end;


figure(110);clf;hold on;
globe([],'earth_1600.png');
arg='100,''yellow'',''fill'',''markeredgecolor'',''black'' ';
for i=1:size(pcapath,3)
    hold on;
    plot3(pcapath(1,:,i),pcapath(2,:,i),pcapath(3,:,i),'y','LineWidth',3);
end
hold on;
plot3(mupath(1,:),mupath(2,:),mupath(3,:),'m','LineWidth',5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% randomly sampling from models
SN = 15; %total sampled data

for i=1:SN
    
%        z = [randn,randn]';
%        diagu = sqrt(diag(Sigmau));
%        sp = z.*diagu;
%        sampled_pca_epmu(:,i)= Up*sp;
%        
%        zw = randn(10,1);
%        diagm = sqrt(diag(Sigmam));
%        spm = zw.*diagm(1:10);
%        sampled_pca_epmw(:,:,i) = reshape(Um(:,1:10)*spm,[T,2])';
    sampled_pca_epmu(:,i) = Um*mvnrnd(mu_new_epmu,Sigmau)';
    sampled_pca_epmw(:,:,i) = reshape(Um*mvnrnd(mu_vnew_epmw,Sigmam)',[T,2])';
end;

%represent the PCA using old coordinate system
for i=1:SN
    sampledoldc_pca_epmu(:,i) = v1*sampled_pca_epmu(1,i) + v2*sampled_pca_epmu(2,i);
    
    for j=1:T
        sampledoldc_pca_epmw(:,j,i) = v1*sampled_pca_epmw(1,j,i) + v2*sampled_pca_epmw(2,j,i);
    end;
end

%reconstruct curves

for i=1:SN
    [sampled_pcap(:,i),sampled_pcaq(:,:,i),sampled_pcapath(:,:,i)] = ExponentialMap(mupath,sampledoldc_pca_epmu(:,i),sampledoldc_pca_epmw(:,:,i),21);
end;






