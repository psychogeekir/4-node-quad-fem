%% Torsion-Prandtl
clc
clear all
close all
%% Load the input file
if ~exist('fem2dinput_TorsionTriangle_Prandtl.mat','file')
    fem2dinput = ConstructInput_Torsion;
else
    load fem2dinput_TorsionTriangle_Prandtl.mat
end

%% Apply the FEM and obtain the stress function
U = fem2dmain(fem2dinput);

%% Some parameters: Torsion Constant, area, polar moment of intertia...
T = 2e3; % real applied torsion
Ex = 2.06e11;
prxy = 0.3;
G = Ex/2/(1+prxy);


A = 0; % Area of cross section
Ixx = 0; % polar moment of inertia = Iyy+Izz, ansys section coordinate
Psi = 0;

for i = 1:fem2dinput.nem
    % begin loop for Gauss Points
    for ngpxi = 1:fem2dinput.ngp
        for ngpeta = 1:fem2dinput.ngp
            % compute the derivatives of shape functions at gauss point
            dpsidxi = fem2dinput.dpsidxi(fem2dinput.gausspt(ngpeta,fem2dinput.ngp));
            dpsideta = fem2dinput.dpsideta(fem2dinput.gausspt(ngpxi,fem2dinput.ngp));
            % compute Jacobian matrix
            J11 = fem2dinput.x(fem2dinput.nod(:,i))'*dpsidxi; 
            J12 = fem2dinput.x(fem2dinput.nod(:,i))'*dpsideta; 
            J21 = fem2dinput.y(fem2dinput.nod(:,i))'*dpsidxi;
            J22 = fem2dinput.y(fem2dinput.nod(:,i))'*dpsideta;
            % compute Jacobian 
            Jac = J11*J22 - J12*J21;
            
            const = fem2dinput.gausswt(ngpxi,fem2dinput.ngp)*fem2dinput.gausswt(ngpeta,fem2dinput.ngp)*Jac;
            
            Ixx = Ixx + ((fem2dinput.x(fem2dinput.nod(:,i))'*...
                fem2dinput.psi(fem2dinput.gausspt(ngpxi,fem2dinput.ngp),fem2dinput.gausspt(ngpeta,fem2dinput.ngp)))^2+...
                (fem2dinput.y(fem2dinput.nod(:,i))'*...
                fem2dinput.psi(fem2dinput.gausspt(ngpxi,fem2dinput.ngp),fem2dinput.gausspt(ngpeta,fem2dinput.ngp)))^2)*const; 
            
            Psi = Psi + U(fem2dinput.nod(:,i))'*...
                fem2dinput.psi(fem2dinput.gausspt(ngpxi,fem2dinput.ngp),fem2dinput.gausspt(ngpeta,fem2dinput.ngp))*const;
            A = A + const;
        end
    end
end
T1 = 2*Psi+2*cell2mat(fem2dinput.HollowArea)'*U(cell2mat(fem2dinput.MasterDOF));
Alpha = T/T1;
J = T1/G;
Uscale = U*Alpha;


%% Plot the stress function
PlotMeshContour(fem2dinput,Uscale,'Torsion_Prandtl.tif')