function Feff = fun_ComputeEquivalentNodalForce(M, C, K, Ue, Uedot, Ueddot, nnodesDRM, InnerDof, OuterDof, t)
%--------------------------------------------------------------------------
% Originally written by Wenyang Zhang
% Email: zwyll@ucla.edu
% Data : 09/12/2018
%
% Refractored by Yiwen Liang
% Email: liangwy@g.ucla.edu
% Data : 12/01/2025
%
% This code computes DRM equivalent nodal forces Feff
%
% INPUT:
%       M, C, K: global mass, damping and stiffness matrices
%
%       Ue, Uedot, Ueddot  : each nnodesDRM*3 by n matrix. displacement, 
%           velocity and acceleration info,
%
%       InnerDof, OuterDof: vectors, node numbers of inner/outer layer of DRM 
%
%       t              : time, a vector.
%
% OUTPUT:
%       Feff           : nnodesDRM*3 by n matrix, equivalent nodel force matrix.
%
%
% Variables:
%       nnodesDRM      : number of nodes in DRM layer (including outter 
%           nodes and inner nodes), also total number of DOFs in the model.
%       n              : length of t.
%       nInnerDof      : length of InnerDof.
%       nOuterDof      : length of OuterDof.
%
%--------------------------------------------------------------------------

n = length(t);

%--------------------------------------------------------------------------
% initialization of the force vector
Mbe = full(M(InnerDof,OuterDof));
Meb = full(M(OuterDof,InnerDof));
Cbe = full(C(InnerDof,OuterDof));
Ceb = full(C(OuterDof,InnerDof));
Kbe = full(K(InnerDof,OuterDof));
Keb = full(K(OuterDof,InnerDof));

% total number of DOFs for the DRM nodes, i.e., DOFs = 3 for each node 
nGlDRMDof = 3*nnodesDRM;

Feff = zeros(nGlDRMDof,n);
linelength = 0;
for i = 1:n
    fprintf(repmat('\b',1,linelength));
    linelength = fprintf('computing equivalent nodal force at time step %d/%d ', i, n);

    Feff(InnerDof,i) = -Mbe*Ueddot(OuterDof,i) - Cbe*Uedot(OuterDof,i) - Kbe*Ue(OuterDof,i);
    Feff(OuterDof,i) = Meb*Ueddot(InnerDof,i) + Ceb*Uedot(InnerDof,i) + Keb*Ue(InnerDof,i);
end 

fprintf('\n');

end