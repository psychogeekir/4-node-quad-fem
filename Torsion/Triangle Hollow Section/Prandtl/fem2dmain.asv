function [U,ReactionOnStructure,R,Q,EBCDOF,NBCDOF] = fem2dmain(fem2dinput)
%% Incorporate the master and slave DOFs


%% Initialize
GlobalK = zeros(fem2dinput.nnodes*fem2dinput.ndof);
GlobalRHS = zeros(fem2dinput.nnodes*fem2dinput.ndof, 1);
GlobalQ = zeros(fem2dinput.nnodes*fem2dinput.ndof, 1);
U = zeros(fem2dinput.nnodes*fem2dinput.ndof,1);

ReactionOnStructure = zeros(fem2dinput.nnodes*fem2dinput.ndof, 1);
R = zeros(fem2dinput.nnodes*fem2dinput.ndof, 1);

elk = zeros(fem2dinput.NodePerElement*fem2dinput.ndof, fem2dinput.NodePerElement*fem2dinput.ndof, fem2dinput.nem);
elrhs = zeros(fem2dinput.NodePerElement*fem2dinput.ndof,fem2dinput.nem);

%% Compute Element K and right-hand-side term
% begin loop for each element
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
            
            % Construct element right-hand-side term, add all gauss points
            elrhs(:,i) = elrhs(:,i) + fem2dinput.f*fem2dinput.psi(fem2dinput.gausspt(ngpxi,fem2dinput.ngp),fem2dinput.gausspt(ngpeta,fem2dinput.ngp))*const;
            
            % Construct element K
            for m = 1:size(elrhs,1)
                for n = 1:size(elrhs,1)
                    elk(m,n,i) = elk(m,n,i) + (fem2dinput.kx*(J11Inv*dpsidxi(m)+J21Inv*dpsideta(m))*(J11Inv*dpsidxi(n)+J21Inv*dpsideta(n))+fem2dinput.ky*(J12Inv*dpsidxi(m)+J22Inv*dpsideta(m))*(J12Inv*dpsidxi(n)+J22Inv*dpsideta(n)))*const;
                end
            end
            
        end
    end
    % end loop for Gaussian Points
end
% end loop for each element

%% Assemble the global K and right-hand-side term
GDOF = unique(fem2dinput.DOFAssociateWithElement(:));
% begin loop for each global DOF
for i = 1:length(GDOF)
    % begin assembling global rhs terms, begin loop for each element
    for k = 1:fem2dinput.nem
        if any(GDOF(i)==fem2dinput.DOFAssociateWithElement(:,k))
            GlobalRHS(i) = GlobalRHS(i) + elrhs(GDOF(i)==fem2dinput.DOFAssociateWithElement(:,k),k);
        end
    end
    % end assembling global rhs terms,  loop for each element
    
    % begin assembling global K, use the symmetry of the K
    for j = i:length(GDOF)
        for n = 1:fem2dinput.nem
            if any(GDOF(i)==fem2dinput.DOFAssociateWithElement(:,n)) && any(GDOF(j)==fem2dinput.DOFAssociateWithElement(:,n))
                GlobalK(i,j) = GlobalK(i,j) + elk(GDOF(i)==fem2dinput.DOFAssociateWithElement(:,n),GDOF(j)==fem2dinput.DOFAssociateWithElement(:,n),n);
            end
        end
    end
    % end assembling global K
end

% Use the symmetry to obtain the complete global K
GlobalK = triu(GlobalK)'+GlobalK-diag(diag(GlobalK));

%% Impose boundary conditions and solve the system
% impose the essential boundary conditions
EBCDOF = false(length(GDOF),1);
EBCDOF(fem2dinput.nodebc) = true;
U(EBCDOF) = fem2dinput.valebc;

% impose the natural boundary conditions
% OldGlobalRHS stores the distributed terms
OldGlobalRHS = GlobalRHS;
NBCDOF = false(length(GDOF),1);
NBCDOF(fem2dinput.nodnbc) = true;
GlobalQ(NBCDOF) = fem2dinput.valnbc;
GlobalRHS(NBCDOF) = GlobalRHS(NBCDOF) + GlobalQ(NBCDOF);

% handle nonzero EBCs
U(~EBCDOF) = 0;
GlobalRHS = GlobalRHS - GlobalK*U;

OldGlobalK = GlobalK;

% reconstruct stiffness matrix using master and slave nodes
if ~isempty(fem2dinput.MasterNode{1,1})
    SlaveDOF = false(length(GDOF),1);
    SlaveDOF(unique(cell2mat(fem2dinput.SlaveDOF))) = true;
    
    for i = 1:length(fem2dinput.MasterDOF)
        GlobalK(fem2dinput.MasterDOF{i,1},:) = GlobalK(fem2dinput.MasterDOF{i,1},:) + sum(GlobalK(fem2dinput.SlaveDOF{i,1},:),1);
        GlobalK(:,fem2dinput.MasterDOF{i,1}) = GlobalK(:,fem2dinput.MasterDOF{i,1}) + sum(GlobalK(:,fem2dinput.SlaveDOF{i,1}),2);      
    end
end

% cross out rows and columns and solve for the global displacement
GlobalK(EBCDOF | SlaveDOF,:) = [];
GlobalK(:,EBCDOF | SlaveDOF) = [];
   
U(~(EBCDOF | SlaveDOF)) = GlobalK\GlobalRHS(~(EBCDOF | SlaveDOF));

% let the results of slave DOFs equal to that of the corresponding master
% DOF
if ~isempty(fem2dinput.MasterNode{1,1})
    for i = 1:length(fem2dinput.MasterDOF)
        U(fem2dinput.SlaveDOF{i,1}) = U(fem2dinput.MasterDOF{i,1});
    end
end

disp('The displacements are')
disp(U)

%% Solve for Secondary Unknowns
Q = OldGlobalK*U - OldGlobalRHS;

% Forces povided by supports
ReactionOnStructure(EBCDOF) = Q(EBCDOF);

% Forces on supports
R(EBCDOF) = -Q(EBCDOF);
R(EBCDOF & NBCDOF) = R(EBCDOF & NBCDOF) + GlobalQ(EBCDOF & NBCDOF);

disp('The reaction forces on supports are')
disp(R(EBCDOF))


end