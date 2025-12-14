clc; clear; close all; tic
%--------------------------------------------------------------------------
% Originally written by Wenyang Zhang
% Email: zwyll@ucla.edu
% Data : 09/12/2018
%
% Refractored by Yiwen Liang
% Email: liangwy@g.ucla.edu
% Data : 11/25/2025
%
% This code generates required DRM information for doing scattered type
% analysis.
% The output is equivalent nodal forces along the DRM interface
%
% before running this code:
% 1) Modify "fun_RayleighDamping"
% 2) Modify "fun_MeshInfor_LSDYNA"
% 3) Modify "fun_NodalResponseAssembly"
%
%--------------------------------------------------------------------------
% this part should be given by the user
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% defining rayleigh damping coefficients
% YOU MUST MODIFY FUNCTION "fun_RayleighDamping"
%--------------------------------------------------------------------------
aC = fun_RayleighDamping; 

%--------------------------------------------------------------------------
% Gauss Quadrature order
QuadOrder = 3;

%--------------------------------------------------------------------------
% reading mesh info from the input file
% YOU MUST MODIFY FUNCTION "fun_MeshInfor_LSDYNA"
%--------------------------------------------------------------------------
[NodeXYCoordinate, ElementConnectivity, NodeData_DRM_inner, NodeData_DRM_outer, ElementData_DRM] = fun_MeshInfor_LSDYNA;

% total number of nodes along the DRM interface
nnodesDRM = length((NodeData_DRM_inner))+length((NodeData_DRM_outer));

% total number of DOFs for the DRM nodes, i.e., DOFs = 3 for each node 
nGlDRMDof = 3*nnodesDRM;

% sorting the DRM nodes' numbers and create a vector for those numbers
DRM_Node_Data = sort([NodeData_DRM_inner;NodeData_DRM_outer]);

% getting the nodes' coordinates for all DRM nodes 
DRM_Node_Coords = find(ismember(NodeXYCoordinate(:,1), DRM_Node_Data));
DRM_Node_Coords = NodeXYCoordinate(DRM_Node_Coords,:);

%--------------------------------------------------------------------------
% computing elemental matrices and assembling
[M, C, K] = fun_Assembly_MCK(DRM_Node_Coords, ElementConnectivity, ElementData_DRM, nnodesDRM, aC, QuadOrder);

disp('done with assembly')

% %-------------------------------------------------------------------------- 
% finding and sorting the inside and outside DRM nodes  
DRM_Inner_Node_Seq = find(ismember(DRM_Node_Coords(:,1), NodeData_DRM_inner));
DRM_Outer_Node_Seq = find(ismember(DRM_Node_Coords(:,1), NodeData_DRM_outer));

% DOFs for inside and outside DRM nodes 
InnerDof = [DRM_Inner_Node_Seq; DRM_Inner_Node_Seq+nnodesDRM; DRM_Inner_Node_Seq+nnodesDRM*2];
OuterDof = [DRM_Outer_Node_Seq; DRM_Outer_Node_Seq+nnodesDRM; DRM_Outer_Node_Seq+nnodesDRM*2];

% %-------------------------------------------------------------------------- 
% compute/load/assemble nodal response for all DRM nodes
% YOU MUST MODIFY FUNCTION "fun_NodalResponseAssembly"
%--------------------------------------------------------------------------
[Ue, Uedot, Ueddot, t] = fun_NodalResponseAssembly(nnodesDRM, DRM_Node_Coords);

disp('done with nodal response assemble')

% %-------------------------------------------------------------------------- 
% compute the equivalent nodal forces for all DRM nodes
Feff = fun_ComputeEquivalentNodalForce(M, C, K, Ue, Uedot, Ueddot, nnodesDRM, InnerDof, OuterDof, t);

disp('done with equivalent nodal force calculation')


%% equivalent force 
% write the files that include the equivalent nodal forces for all DRM
%   nodes in 3 directions (Fx, Fy, Fz). you can save or load required
%   varibales from EquF.mat from previous runs.

tic
load EquF.mat
 
disp('exporting equivalent force to LS-DYNA keyword file...')

fun_WriteNodalForcestoLSDYNAinput(nnodesDRM,DRM_Node_Data,t,Feff,10000);

%   save('EquF.mat', 'nnodesDRM', 'DRM_Node_Data', 't', 'Feff', '-v7.3');

disp('program finished')
toc

