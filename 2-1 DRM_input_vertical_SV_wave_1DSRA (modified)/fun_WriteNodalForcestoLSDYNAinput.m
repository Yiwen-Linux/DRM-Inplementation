function fun_WriteNodalForcestoLSDYNAinput(nnodesDRM,DRM_Node_Data,t,Feff)
%--------------------------------------------------------------------------
% Originally written by Wenyang Zhang
% Email: zwyll@ucla.edu
% Data : 09/12/2018
%
% Modified by Yiwen
% Data : 11/25/2025
%
% This code writes DRM equivalent nodal forces to .txt file
%
% INPUT:
%       nnodesDRM      : number of nodes in DRM layer (including outter 
%           nodes and inner nodes.
%       DRM_Node_Data  : node numbers of DRM layer, a nnodesDRM by 1
%           matrix.
%       t              : time, a vector.
%       Feff           : equivalent nodal forces at each node, a nnodesDRM
%           by n matrix, n is the length of t.
%
% OUTPUT:
%       (NONE)
%
%
%--------------------------------------------------------------------------


mkdir('Files need to be imported into LS-DYNA');
n = length(t);

for i = 1:nnodesDRM
    OutputName = sprintf('Files need to be imported into LS-DYNA/EqNodalForceInfo%d_x.txt',i);
    fid = fopen(OutputName,'wt');
    fprintf(fid,'*DEFINE_CURVE\n');
    fprintf(fid,'%d\n',i);
    for j = 1:n
        fprintf(fid,'%.4f,   %.6f\n',t(j),Feff(i,j));
    end
    fclose(fid);
    
    OutputName = sprintf('Files need to be imported into LS-DYNA/EqNodalForceInfo%d_y.txt',i);
    fid = fopen(OutputName,'wt');
    fprintf(fid,'*DEFINE_CURVE\n');
    fprintf(fid,'%d\n',i+nnodesDRM);
    for j = 1:n
        fprintf(fid,'%.4f,   %.6f\n',t(j),Feff(i+nnodesDRM,j));
    end
    fclose(fid);
    
    OutputName = sprintf('Files need to be imported into LS-DYNA/EqNodalForceInfo%d_z.txt',i);
    fid = fopen(OutputName,'wt');
    fprintf(fid,'*DEFINE_CURVE\n');
    fprintf(fid,'%d\n',i+nnodesDRM*2);
    for j = 1:n
        fprintf(fid,'%.4f,   %.6f\n',t(j),Feff(i+nnodesDRM*2,j));
    end
    fclose(fid);      
end

% writing the LS-DYNA keywor files to include the files for equivalent
% nodal forces curves
fid = fopen('Files need to be imported into LS-DYNA/IncludeInfo.txt','wt');
for i = 1:nnodesDRM 
    fprintf(fid,'*INCLUDE\n');
    fprintf(fid,' EqNodalForceInfo%d_x.txt\n',i);
    fprintf(fid,'*INCLUDE\n');
    fprintf(fid,' EqNodalForceInfo%d_y.txt\n',i);
    fprintf(fid,'*INCLUDE\n');
    fprintf(fid,' EqNodalForceInfo%d_z.txt\n',i);
end

% writing the LS-DYNA files to define nodal forces as *LOAD_NODE_POINT
fid = fopen('Files need to be imported into LS-DYNA/LoadInfo.txt','wt');
fprintf(fid,'*LOAD_NODE_POINT\n');
for i = 1:nnodesDRM  
    fprintf(fid,'%d,1,%d,1.0,0,0,0,0\n',DRM_Node_Data(i),i);
    fprintf(fid,'%d,2,%d,1.0,0,0,0,0\n',DRM_Node_Data(i),i+nnodesDRM);
    fprintf(fid,'%d,3,%d,1.0,0,0,0,0\n',DRM_Node_Data(i),i+nnodesDRM*2);
end
fclose(fid);

end