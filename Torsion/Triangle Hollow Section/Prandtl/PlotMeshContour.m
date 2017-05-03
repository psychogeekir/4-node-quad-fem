function PlotMeshContour(fem2dinput,U,PicName)

%% Plot

figure
for iElement = 1:fem2dinput.nem
    xNode = fem2dinput.x(fem2dinput.nod(:,iElement));
    yNode = fem2dinput.y(fem2dinput.nod(:,iElement));
    % fill3(xNode,yNode,ones(length(xNode),1),'w','LineWidth',0.25)
    fill(xNode,yNode,U(fem2dinput.nod(:,iElement)),'LineWidth',0.25)
    %alpha 0
    hold on
end

% [X,Y] = meshgrid(fem2dinput.x,fem2dinput.y);
% [Xq,Yq,Zq] = griddata(fem2dinput.x,fem2dinput.y,U,X,Y);
% [C,t] = contour(Xq,Yq,Zq);
% t.LineWidth = 0.25;
% t.LineColor = 'red';
% t.LevelStep = 3;
% t.Fill = 'on';
% clabel(C,t)
% view(2)
hold off

axis([0 fem2dinput.TotalDomain(1) 0 fem2dinput.TotalDomain(2)])
axis equal
axis tight
colorbar
xlabel('\it x')
ylabel('\it y')
title(fem2dinput.PlotTitle)
h=gca;
% h.XTick = unique(fem2dinput.x);
% h.YTick = unique(fem2dinput.y);
h.FontName='Times New Roman';
h.FontSize=11; 
set(gcf,'Position',[400 400 500 400],'Color','w')
export_fig(PicName,'-r600','-opengl')

% figure
% surf(X,Y,Z,'EdgeColor','None','facecolor','interp')
% colorbar
% xlabel('\it x')
% ylabel('\it y')
% zlabel('\it U')
% title(fem2dinput.PlotTitle)
% h=gca;
% h.XTick = unique(fem2dinput.x);
% h.YTick = unique(fem2dinput.y);
% grid on
% grid minor
% h.FontName='Times New Roman';
% h.FontSize=11; 
% set(gcf,'Position',[400 400 450 350],'Color','w')
% export_fig([fem2dinput.PlotTitle,'_3d.tif'],'-r600','-opengl')

end