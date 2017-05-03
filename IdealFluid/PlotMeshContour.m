function PlotMeshContour(fem2dinput,U,PicName)

%% Plot
% Define Shape Functions for the element parallel with the global coordinate
N = @(x,y,xNode,yNode,Length) 1/4*(1+2*sign(xNode).*x/Length(1)).*(1+2*sign(yNode).*y/Length(2));

[X,Y] = meshgrid(0:0.01:4,0:0.01:2);
[dimi,dimj] = size(X);
Z = NaN(size(X));
for i = 1:dimi
    for j = 1:dimj    
        if ~(X(i,j)>=3 && X(i,j)<=4 && Y(i,j)>=0 && Y(i,j)<=1)
            [NodeNumber,ElementNumber] = GetCurrentElement(X(i,j),Y(i,j),fem2dinput.nod,[fem2dinput.x,fem2dinput.y]);
            % Obtain the coordinates of the adjacent nodes and translate
            xNode = fem2dinput.x(fem2dinput.nod(:,ElementNumber));
            yNode = fem2dinput.y(fem2dinput.nod(:,ElementNumber));
            xmean = mean(xNode);
            ymean = mean(yNode);
            xNode = xNode - xmean;
            yNode = yNode - ymean;
            % Compute the displacement using interpolation
            Z(i,j) = N(X(i,j) - xmean,Y(i,j) - ymean,xNode,yNode,fem2dinput.ElementSize(ElementNumber,:))'*U(fem2dinput.DOFAssociateWithElement(:,ElementNumber));   
        end
    end
end

figure
for iElement = 1:fem2dinput.nem
    xNode = fem2dinput.x(fem2dinput.nod(:,iElement));
    yNode = fem2dinput.y(fem2dinput.nod(:,iElement));
    fill3(xNode,yNode,ones(length(xNode),1),'w','LineWidth',0.25)
    alpha 0
    hold on
end

[C,t] = contour(X,Y,Z);
t.LineWidth = 2;
t.LineColor = 'red';
t.LevelStep = 0.25;
t.Fill = 'on';
clabel(C,t)
view(2)
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