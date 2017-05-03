%% Ideal Fluid Potential
% clc
% clear all
% close all
% 
% %%
% if ~exist('fem2dinput_Potential_0.2x0.2.mat','file')
%     fem2dinput = ConstructInput_Potential;
% else
%     load fem2dinput_Potential_0.2x0.2.mat
% end
% 
% %%
% U = fem2dmain(fem2dinput);
% 
% %%
% PlotMeshContour(fem2dinput,U,'Ideal Fluid Potential.tif')

%% Ideal Fluid Streamline
clc
clear all
close all
%%
if ~exist('fem2dinput_Streamline_0.2x0.2.mat','file')
    fem2dinput = ConstructInput_Streamline;
else
    load fem2dinput_Streamline_0.2x0.2.mat
end

%%
U = fem2dmain(fem2dinput);

%%
PlotMeshContour(fem2dinput,U,'Ideal Fluid Streamline.tif')