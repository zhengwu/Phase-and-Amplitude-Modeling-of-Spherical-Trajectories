function PLOT(X,arg,color,linewidth)
% [sphereX,sphereY,sphereZ] =sphere;
r=.96;
r2=.98;
LW = linewidth;

if ~exist('arg','var')
    arg='100,''fill'',''markeredgecolor'',''black'' ';
end
if arg==1
    plot3(X(1,:),X(2,:),X(3,:),color,'linewidth',3)
    axis equal;set(gca,'linewidth',LW,'fontsize',18)
else
eval(['scatter3(X(1,:),X(2,:),X(3,:),' arg ')'])
axis equal;set(gca,'linewidth',2,'fontsize',18)
end
end
