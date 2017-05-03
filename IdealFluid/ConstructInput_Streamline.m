function fem2dinput = ContructInput_Streamline

clear variables

%% 0.2x0.2
fem2dinput.PlotTitle = 'Ideal Fluid -Streamline -0.2\times0.2 Element';

fem2dinput.TotalDomain = [4, 2]; % The rectangle wrapping the whole domain

fem2dinput.nem = (3/0.2)*(2/0.2)+(1/0.2)*(1/0.2);
fem2dinput.ElementSize = 0.2*ones(fem2dinput.nem,2);
fem2dinput.nnodes = (3/0.2+1)*(2/0.2+1)+(1/0.2+1)*(1/0.2);
fem2dinput.NodePerElement = 4;
fem2dinput.ndof = 1;

fem2dinput.kx = 1;
fem2dinput.ky = 1;
fem2dinput.f = 0;

%%
fem2dinput.x = zeros(fem2dinput.nnodes,1);
fem2dinput.y = zeros(fem2dinput.nnodes,1);
Esize = fem2dinput.ElementSize(1,1);
for i = 1:(3/0.2+1)
    fem2dinput.x(i:(3/0.2+1):(i+2/0.2*(3/0.2+1))) = 0 + 0.2*(i-1);
    fem2dinput.y(i:(3/0.2+1):(i+2/0.2*(3/0.2+1))) = 0:0.2:2;
end
for j = 1:(1/0.2)
    fem2dinput.x((i+2/0.2*(3/0.2+1)+1+(j-1)*(1/0.2+1)):(i+2/0.2*(3/0.2+1)+1+1/0.2+(j-1)*(1/0.2+1))) = 3+j*0.2;
    fem2dinput.y((i+2/0.2*(3/0.2+1)+1+(j-1)*(1/0.2+1)):(i+2/0.2*(3/0.2+1)+1+1/0.2+(j-1)*(1/0.2+1))) = 1:0.2:2;
end

%%
fem2dinput.nod = zeros(fem2dinput.NodePerElement, fem2dinput.nem);
centerx = zeros(fem2dinput.nem,1);
centery = zeros(fem2dinput.nem,1);
for i = 1:(3/0.2)
    centerx(i:(3/0.2):(i+(2/0.2-1)*(3/0.2))) = 0.2/2 + 0.2*(i-1);
    centery(i:(3/0.2):(i+(2/0.2-1)*(3/0.2))) = (0+0.2/2):0.2:(2-0.2/2);
end
for j = 1:(1/0.2)
    centerx((i+(2/0.2-1)*(3/0.2)+1+(j-1)*(1/0.2)):(i+(2/0.2-1)*(3/0.2)+1+1/0.2-1+(j-1)*(1/0.2))) = 3+0.2/2+(j-1)*0.2;
    centery((i+(2/0.2-1)*(3/0.2)+1+(j-1)*(1/0.2)):(i+(2/0.2-1)*(3/0.2)+1+1/0.2-1+(j-1)*(1/0.2))) = (1+0.2/2):0.2:(2-0.2/2);
end

for i = 1:fem2dinput.nem
   fem2dinput.nod(1,i) = find(abs(centerx(i)-0.2/2-fem2dinput.x)<=1e-6 & abs(centery(i)-0.2/2-fem2dinput.y)<=1e-7);
   fem2dinput.nod(2,i) = find(abs(centerx(i)+0.2/2-fem2dinput.x)<=1e-6 & abs(centery(i)-0.2/2-fem2dinput.y)<=1e-7);
   fem2dinput.nod(3,i) = find(abs(centerx(i)+0.2/2-fem2dinput.x)<=1e-6 & abs(centery(i)+0.2/2-fem2dinput.y)<=1e-7);
   fem2dinput.nod(4,i) = find(abs(centerx(i)-0.2/2-fem2dinput.x)<=1e-6 & abs(centery(i)+0.2/2-fem2dinput.y)<=1e-7);
end

fem2dinput.DOFAssociateWithElement = fem2dinput.nod;


%%
nodebc{1,1} = find(abs(0-fem2dinput.x)<=1e-6);
valebc{1,1} = fem2dinput.y(nodebc{1,1});

nodebc{2,1} = find(abs(0-fem2dinput.y)<=1e-6);
valebc{2,1} = zeros(length(nodebc{2,1}),1);

nodebc{3,1} = find(abs(3-fem2dinput.x)<=1e-6 & fem2dinput.y>=0 & fem2dinput.y<=1);
valebc{3,1} = zeros(length(nodebc{3,1}),1);

nodebc{4,1} = find(abs(1-fem2dinput.y)<=1e-6 & fem2dinput.x>=3 & fem2dinput.x<=4);
valebc{4,1} = zeros(length(nodebc{4,1}),1);

nodebc{5,1} = find(abs(2-fem2dinput.y)<=1e-6);
valebc{5,1} = ones(length(nodebc{5,1}),1)*fem2dinput.y(abs(0-fem2dinput.x)<=1e-6 & abs(2-fem2dinput.y)<=1e-6);


fem2dinput.nodebc = unique(cell2mat(nodebc));
fem2dinput.valebc = zeros(length(fem2dinput.nodebc),1);
for i = 1:length(fem2dinput.nodebc)
    for j = 1:size(nodebc,1)
        if ~isempty(valebc{j,1}(nodebc{j,1} == fem2dinput.nodebc(i))) 
            fem2dinput.valebc(i) = valebc{j,1}(nodebc{j,1} == fem2dinput.nodebc(i));
        end
    end
end
fem2dinput.nebc = length(fem2dinput.nodebc);

%%
nodnbc{1,1} = find(abs(4-fem2dinput.x)<=1e-6);
valnbc{1,1} = zeros(length(nodnbc{1,1}),1);

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
fem2dinput.ngp = 2;

fem2dinput.gausspt = zeros(4,4);
fem2dinput.gausswt = zeros(4,4);
for i = 1:size(fem2dinput.gausswt,1)
    [fem2dinput.gausspt(1:i,i), fem2dinput.gausswt(1:i,i)] = md_gauss(i,1);
end

fem2dinput.psi = @(xi,eta) [1/4*(1-xi)*(1-eta);...
                            1/4*(1+xi)*(1-eta);...
                            1/4*(1+xi)*(1+eta);...
                            1/4*(1-xi)*(1+eta)];

fem2dinput.dpsidxi = @(eta) [-1/4*(1-eta);1/4*(1-eta);1/4*(1+eta);-1/4*(1+eta)];

fem2dinput.dpsideta = @(xi) [-1/4*(1-xi);-1/4*(1+xi);1/4*(1+xi);1/4*(1-xi)];

save('fem2dinput_Streamline_0.2x0.2.mat','fem2dinput')

end
