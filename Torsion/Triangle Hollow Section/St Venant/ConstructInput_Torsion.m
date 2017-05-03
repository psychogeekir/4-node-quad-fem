function fem2dinput = ConstructInput_Torsion

clear variables

%% 
fem2dinput.PlotTitle = 'Torsion -St Venant';

fem2dinput.TotalDomain = [10, 6]; % The rectangle wrapping the whole domain

fem2dinput.ElementSize = [];

fem2dinput.NodePerElement = 4;
fem2dinput.ndof = 1;

fem2dinput.kx = 1;
fem2dinput.ky = 1;
fem2dinput.f = 0;

%%
temp = importdata('NodeCoordinate.dat',' ',4);
fem2dinput.nnodes = size(temp.data,1);
fem2dinput.x = temp.data(:,2);
fem2dinput.y = temp.data(:,3);

% Note: in Ansys, the node associated with the element is in clockwise (negative Jacobian) or countclockwise direction (positive Jacobian), 
% need to transform it to countclockwise direction (positive Jacobian), consistent with the order of shape function 
temp = importdata('NodeAssociatedElement.dat',' ',5);
fem2dinput.nem = size(temp.data,1);
fem2dinput.nod = zeros(fem2dinput.NodePerElement,fem2dinput.nem);
for i = 1:fem2dinput.nem
    fem2dinput.nod(:,i) = ReorderNumberingSequence(temp.data(i,end-3:end)',[fem2dinput.x(temp.data(i,end-3:end)'),fem2dinput.y(temp.data(i,end-3:end)')]);
end

clear temp

fem2dinput.DOFAssociateWithElement = fem2dinput.nod;



%%
fem2dinput.ngp = 2;

fem2dinput.gausspt = zeros(4,4);
fem2dinput.gausswt = zeros(4,4);
for i = 1:size(fem2dinput.gausswt,1)
    [fem2dinput.gausspt(1:i,i), fem2dinput.gausswt(1:i,i)] = md_gauss(i,1);
end

% anticlockwise
fem2dinput.psi = @(xi,eta) [1/4*(1-xi)*(1-eta);...
                            1/4*(1+xi)*(1-eta);...
                            1/4*(1+xi)*(1+eta);...
                            1/4*(1-xi)*(1+eta)];

fem2dinput.dpsidxi = @(eta) [-1/4*(1-eta);1/4*(1-eta);1/4*(1+eta);-1/4*(1+eta)];

fem2dinput.dpsideta = @(xi) [-1/4*(1-xi);-1/4*(1+xi);1/4*(1+xi);1/4*(1-xi)];


%%
% Just to remove the singularity
fem2dinput.nodebc = [1];
fem2dinput.valebc = [0];
fem2dinput.nebc = length(fem2dinput.nodebc);

%%
psi = {@(ksi) (1-ksi)/2; @(ksi) (1+ksi)/2};
% qn = y*nx-x*ny
nxy = [1,0; cos(120/180*pi), sin(120/180*pi); cos(240/180*pi), sin(240/180*pi);...
      -1,0; -cos(120/180*pi), -sin(120/180*pi); -cos(240/180*pi), -sin(240/180*pi)];

qn = cell(length(nxy),1); 
for i = 1:length(nxy)
    qn{i,1} = @(x,y) [y,-x]*nxy(i,:)';
end


nodnbc{1,1} = find(abs(1/3*sqrt(3)/2*10-fem2dinput.x)<=1e-6);
valnbc{1,1} = zeros(length(nodnbc{1,1}),1);
valnbc{1,1} = CalLineNBC(nodnbc{1,1},valnbc{1,1},[fem2dinput.x(nodnbc{1,1}),fem2dinput.y(nodnbc{1,1})],'y',psi,qn{1},fem2dinput.gausspt,fem2dinput.gausswt);

nodnbc{2,1} = find(abs(tan(pi/6)*(fem2dinput.x+2/3*sqrt(3)/2*10)-fem2dinput.y)<=1e-6);
valnbc{2,1} = zeros(length(nodnbc{2,1}),1);
valnbc{2,1} = CalLineNBC(nodnbc{2,1},valnbc{2,1},[fem2dinput.x(nodnbc{2,1}),fem2dinput.y(nodnbc{2,1})],'x',psi,qn{2},fem2dinput.gausspt,fem2dinput.gausswt);

nodnbc{3,1} = find(abs(-tan(pi/6)*(fem2dinput.x+2/3*sqrt(3)/2*10)-fem2dinput.y)<=1e-6);
valnbc{3,1} = zeros(length(nodnbc{3,1}),1);
valnbc{3,1} = CalLineNBC(nodnbc{3,1},valnbc{3,1},[fem2dinput.x(nodnbc{3,1}),fem2dinput.y(nodnbc{3,1})],'x',psi,qn{3},fem2dinput.gausspt,fem2dinput.gausswt);

nodnbc{4,1} = find(abs(1/3*sqrt(3)/2*10-0.5-fem2dinput.x)<=1e-6 & fem2dinput.y>=-10/2+sqrt(3)*0.5-0.1 & fem2dinput.y<=10/2-sqrt(3)*0.5+0.1);
valnbc{4,1} = zeros(length(nodnbc{4,1}),1);
valnbc{4,1} = CalLineNBC(nodnbc{4,1},valnbc{4,1},[fem2dinput.x(nodnbc{4,1}),fem2dinput.y(nodnbc{4,1})],'y',psi,qn{4},fem2dinput.gausspt,fem2dinput.gausswt);

nodnbc{5,1} = find(abs(tan(pi/6)*(fem2dinput.x+2/3*sqrt(3)/2*10-1)-fem2dinput.y)<=1e-6 & fem2dinput.y>=0-0.1 & fem2dinput.y<=10/2-sqrt(3)*0.5+0.1);
valnbc{5,1} = zeros(length(nodnbc{5,1}),1);
valnbc{5,1} = CalLineNBC(nodnbc{5,1},valnbc{5,1},[fem2dinput.x(nodnbc{5,1}),fem2dinput.y(nodnbc{5,1})],'x',psi,qn{5},fem2dinput.gausspt,fem2dinput.gausswt);

nodnbc{6,1} = find(abs(-tan(pi/6)*(fem2dinput.x+2/3*sqrt(3)/2*10-1)-fem2dinput.y)<=1e-6 & fem2dinput.y<=0+0.1 & fem2dinput.y>=-10/2+sqrt(3)*0.5-0.1);
valnbc{6,1} = zeros(length(nodnbc{6,1}),1);
valnbc{6,1} = CalLineNBC(nodnbc{6,1},valnbc{6,1},[fem2dinput.x(nodnbc{6,1}),fem2dinput.y(nodnbc{6,1})],'x',psi,qn{6},fem2dinput.gausspt,fem2dinput.gausswt);


fem2dinput.nodnbc = unique(cell2mat(nodnbc));
fem2dinput.valnbc = zeros(length(fem2dinput.nodnbc),1);
for i = 1:length(fem2dinput.nodnbc)
    for j = 1:size(nodnbc,1)
        if ~isempty(valnbc{j,1}(nodnbc{j,1} == fem2dinput.nodnbc(i))) 
            fem2dinput.valnbc(i) = fem2dinput.valnbc(i) + valnbc{j,1}(nodnbc{j,1} == fem2dinput.nodnbc(i));
        end
    end
end
fem2dinput.nnbc = length(fem2dinput.nodnbc);

%%
save('fem2dinput_TorsionTriangle_StVenant.mat','fem2dinput')

end

function CCWNodeNo = ReorderNumberingSequence(NodeNo,Coordinate)
% Use cross product to check if the sequence is in counter-clockwise direction
    
    if size(Coordinate,2)==2
        Coordinate = [Coordinate, zeros(size(Coordinate,1),1)];
    end
    
    AdjacentVector = diff([Coordinate(end,:); Coordinate; Coordinate(1,:)]);
    tmp = cross(AdjacentVector(1:end-1,:)',AdjacentVector(2:end,:)')';
    if all(sign(tmp(:,3)) == 1)
        CCWNodeNo = NodeNo;
    elseif all(sign(tmp(:,3)) == -1)
        CCWNodeNo = flipud(NodeNo);
    else
        s = randperm(length(NodeNo));
        CCWNodeNo = ReorderNumberingSequence(NodeNo(s,1),Coordinate(s,:));
    end
end

function Q = Cal1DGQ(qn,psi,xpq,ypq,gausspt,gausswt)
% xp<=xq
% The line is represented by parametric form (vector) 
% xx = (1-t)*xp + t*xq, yy = (1-t)*yp + t*yq, t=0..1
t = @(ksi) (ksi*1+0+1)/2;
xx = @(t) (1-t)*xpq(1) + t*xpq(2);
yy = @(t) (1-t)*ypq(1) + t*ypq(2);
Q = qn(xx(t(gausspt(1,2))),yy(t(gausspt(1,2))))*psi(gausspt(1,2))*sqrt((xpq(2)-xpq(1))^2+(ypq(2)-ypq(1))^2)*gausswt(1,2)*1/2 ...
    + qn(xx(t(gausspt(2,2))),yy(t(gausspt(2,2))))*psi(gausspt(2,2))*sqrt((xpq(2)-xpq(1))^2+(ypq(2)-ypq(1))^2)*gausswt(2,2)*1/2;

end

function valnbc = CalLineNBC(nodnbc,valnbc,coordinate,sortcriterion,psi,qn,gausspt,gausswt)
if sortcriterion=='x';
    tt = sortrows([nodnbc,coordinate],2);
elseif sortcriterion=='y'
    tt = sortrows([nodnbc,coordinate],3);
end

for i = 1:length(nodnbc)
    if i == 1
        valnbc(nodnbc == tt(i,1)) = Cal1DGQ(qn,psi{1},[tt(i,2),tt(i+1,2)],[tt(i,3),tt(i+1,3)],gausspt,gausswt);
    elseif i==length(nodnbc)
        valnbc(nodnbc == tt(i,1)) = Cal1DGQ(qn,psi{2},[tt(i-1,2),tt(i,2)],[tt(i-1,3),tt(i,3)],gausspt,gausswt);
    else
        valnbc(nodnbc == tt(i,1)) = Cal1DGQ(qn,psi{1},[tt(i,2),tt(i+1,2)],[tt(i,3),tt(i+1,3)],gausspt,gausswt)+...
            Cal1DGQ(qn,psi{2},[tt(i-1,2),tt(i,2)],[tt(i-1,3),tt(i,3)],gausspt,gausswt);
    end            
end

end
