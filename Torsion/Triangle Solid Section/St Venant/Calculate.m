%% Torsion-St Venant
clc
clear all
close all

%% Load the input file
if ~exist('fem2dinput_TorsionTriangle_StVenant.mat','file')
    fem2dinput = ConstructInput_Torsion;
else
    load fem2dinput_TorsionTriangle_StVenant.mat
end

%% Apply the FEM and obtain the warping function
U = fem2dmain(fem2dinput);

%% Some parameters: Torsion Constant, area, polar moment of intertia...
T = 2e3;

Ex = 2.06e11;
prxy = 0.3;
G = Ex/2/(1+prxy);

A = 0; % area
J = 0; % torsion constant
Ixx = 0; % polar moment of inertia = Iyy+Izz, ansys section coordinate

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
            % compute the inverse of Jacobian matrix
            J11Inv = J22/Jac;
            J22Inv = J11/Jac;
            J12Inv = -J12/Jac;
            J21Inv = -J21/Jac;
            
            const = fem2dinput.gausswt(ngpxi,fem2dinput.ngp)*fem2dinput.gausswt(ngpeta,fem2dinput.ngp)*Jac;

            Ixx = Ixx + ((fem2dinput.x(fem2dinput.nod(:,i))'*...
                fem2dinput.psi(fem2dinput.gausspt(ngpxi,fem2dinput.ngp),fem2dinput.gausspt(ngpeta,fem2dinput.ngp)))^2+...
                (fem2dinput.y(fem2dinput.nod(:,i))'*...
                fem2dinput.psi(fem2dinput.gausspt(ngpxi,fem2dinput.ngp),fem2dinput.gausspt(ngpeta,fem2dinput.ngp)))^2)*const;
            
            J = J +...
                ((fem2dinput.x(fem2dinput.nod(:,i))'*...
                fem2dinput.psi(fem2dinput.gausspt(ngpxi,fem2dinput.ngp),fem2dinput.gausspt(ngpeta,fem2dinput.ngp)))*...
                (U(fem2dinput.nod(:,i))'*dpsidxi*J12Inv+U(fem2dinput.nod(:,i))'*dpsideta*J22Inv)-...
                (fem2dinput.y(fem2dinput.nod(:,i))'*...
                fem2dinput.psi(fem2dinput.gausspt(ngpxi,fem2dinput.ngp),fem2dinput.gausspt(ngpeta,fem2dinput.ngp)))*...
                (U(fem2dinput.nod(:,i))'*dpsidxi*J11Inv+U(fem2dinput.nod(:,i))'*dpsideta*J21Inv))*const;
            
            A = A + const;
        end
    end
end

J = Ixx+J;
Alpha = T/G/J;
%% Plot the warping function
PlotMeshContour(fem2dinput,U,'Torsion_StVenant.tif')

