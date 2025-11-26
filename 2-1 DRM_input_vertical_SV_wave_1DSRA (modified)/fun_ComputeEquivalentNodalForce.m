function Feff = fun_ComputeEquivalentNodalForce(Mbe, Meb, Cbe, Ceb, Kbe, Keb, Ue, Uedot, Ueddot, InnerDof, OuterDof, t)
%--------------------------------------------------------------------------
% Originally written by Wenyang Zhang
% Email: zwyll@ucla.edu
% Data : 09/12/2018
%
% Refractored by Yiwen
% Data : 11/25/2025
%
% This code computes DRM equivalent nodal forces Feff
%
% INPUT:
%       Mbe, Cbe, Kbe: nInnerDof by n OuterDof matrices with
%           subscripts*, partioning of mass, damping and stiffness matrices
%           in b and e crossing area, 
%       Meb, Ceb, Keb: OuterDof by n nInnerDof matrices with
%           subscripts*, partioning of mass, damping and stiffness matrices
%           in e and b crossing area, 
%
%           * the subscripts b and e refer to the nodes along the inside
%             and outside boundary of the one layer of elements,
%             respectively.
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

Feff = zeros(nGlDRMDof,n);
for i = 1:n
    Feff(InnerDof,i) = -Mbe*Ueddot(OuterDof,i) - Cbe*Uedot(OuterDof,i) - Kbe*Ue(OuterDof,i);
    Feff(OuterDof,i) = Meb*Ueddot(InnerDof,i) + Ceb*Uedot(InnerDof,i) + Keb*Ue(InnerDof,i);
end 

end