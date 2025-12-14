function fun_WriteNodalForcestoLSDYNAinput(nnodesDRM,DRM_Node_Data,t,Feff,countn_forcecurve)
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
%       countn_forcecurve: scalar, starting point of force curve
%
% OUTPUT:
%       (NONE)
%
%
%--------------------------------------------------------------------------


mkdir('Files need to be imported into LS-DYNA');
n = length(t);


% writing the LS-DYNA keywor files to include the files for equivalent
% nodal forces curves
disp('writing LS-DYNA include files')

fid = fopen('Files need to be imported into LS-DYNA/IncludeInfo.txt','wt');
for i = 1:nnodesDRM 
    fprintf(fid,'*INCLUDE\n');
    fprintf(fid,'Files need to be imported into LS-DYNA\\EqNodalForceInfo%d_x.txt\n',i);
    fprintf(fid,'*INCLUDE\n');
    fprintf(fid,'Files need to be imported into LS-DYNA\\EqNodalForceInfo%d_y.txt\n',i);
    fprintf(fid,'*INCLUDE\n');
    fprintf(fid,'Files need to be imported into LS-DYNA\\EqNodalForceInfo%d_z.txt\n',i);
end

disp('writing LS-DYNA include files completed');

fclose(fid);

% writing the LS-DYNA files to define nodal forces as *LOAD_NODE_POINT
disp('writing LS-DYNA nodal force definition')

fid = fopen('Files need to be imported into LS-DYNA/LoadInfo.txt','wt');
fprintf(fid,'*LOAD_NODE_POINT\n');
for i = 1:nnodesDRM  
    fprintf(fid,'%d,1,%d,1.0,0,0,0,0\n',DRM_Node_Data(i),countn_forcecurve+i);
    fprintf(fid,'%d,2,%d,1.0,0,0,0,0\n',DRM_Node_Data(i),countn_forcecurve+i+nnodesDRM);
    fprintf(fid,'%d,3,%d,1.0,0,0,0,0\n',DRM_Node_Data(i),countn_forcecurve+i+nnodesDRM*2);
end

disp('LS-DYNA nodal force definition completed')

fclose(fid);

% writing the equivalent nodal force time seires as *DEFINE_CURVE
linelength = 0;
for i = 1:nnodesDRM
    fprintf(repmat('\b',1,linelength));
    linelength = fprintf('writing equivalent nodal force %d/%d', i, nnodesDRM);

    OutputName = sprintf('Files need to be imported into LS-DYNA/EqNodalForceInfo%d_x.txt',i);
    fid = fopen(OutputName,'wt');
    fprintf(fid,'*DEFINE_CURVE\n');
    fprintf(fid,'%d\n',countn_forcecurve+i);
    for j = 1:n
        fprintf(fid,'%.4f,   %.6f\n',t(j),Feff(i,j));
    end
    fclose(fid);
    
    OutputName = sprintf('Files need to be imported into LS-DYNA/EqNodalForceInfo%d_y.txt',i);
    fid = fopen(OutputName,'wt');
    fprintf(fid,'*DEFINE_CURVE\n');
    fprintf(fid,'%d\n',countn_forcecurve+i+nnodesDRM);
    for j = 1:n
        fprintf(fid,'%.4f,   %.6f\n',t(j),Feff(i+nnodesDRM,j));
    end
    fclose(fid);
    
    OutputName = sprintf('Files need to be imported into LS-DYNA/EqNodalForceInfo%d_z.txt',i);
    fid = fopen(OutputName,'wt');
    fprintf(fid,'*DEFINE_CURVE\n');
    fprintf(fid,'%d\n',countn_forcecurve+i+nnodesDRM*2);
    for j = 1:n
        fprintf(fid,'%.4f,   %.6f\n',t(j),Feff(i+nnodesDRM*2,j));
    end
    fclose(fid);      
end

fprintf('\n');
disp('writing equivalent nodal force files completed')
end