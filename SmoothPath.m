function smoothed_path=SmoothPath(path,WindSize,Sig)
% Gaussian mean filter
if(mod(WindSize,2)~=1)
    WindSize  = WindSize + 1;
end
lN = (WindSize -1)/2;

filter = fspecial('gaussian',[1,WindSize],Sig);


[n,T]=size(path);
for i = 1:T
    %smoothing
    Tmp_mean_point = zeros(n,1);
    wid = 1;
    tweight = 0;
    for j=-lN:lN
        idx = i+j;
        if(idx>0 && idx<T+1)
            Tmp_mean_point = Tmp_mean_point + path(:,idx)*filter(wid);
            tweight = tweight+filter(wid);
        end;
        wid = wid + 1;
    end;
    tmp_p =  Tmp_mean_point/tweight;
    smoothed_path(:,i) = tmp_p/norm(tmp_p,'fro');
end;
