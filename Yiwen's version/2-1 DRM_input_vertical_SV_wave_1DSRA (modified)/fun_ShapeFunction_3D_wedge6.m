function [N, B] = fun_ShapeFunction_3D_wedge6(xi, eta, mu)
% fun_ShapeFunction_3D_wedge6
% 6-node linear wedge element (triangle x line)
%
% Input:
%   xi, eta, mu : natural coordinates
%       - triangle coordinates: xi >= 0, eta >= 0, xi + eta <= 1
%       - thickness coordinate: mu in [-1, 1]
%
% Output:
%   N : 6x1 shape function vector
%   B : 3x6 matrix, each column = [dNi/dxi; dNi/deta; dNi/dmu]

    % barycentric coordinates on the triangle
    L1 = 1 - xi - eta;
    L2 = xi;
    L3 = eta;

    % vertical interpolation factors
    a = (1 - mu)/2;   % bottom (mu = -1 -> a=1, b=0)
    b = (1 + mu)/2;   % top    (mu =  1 -> a=0, b=1)

    % shape functions (6 by 1)
    N1 = L1 * a;
    N2 = L2 * a;
    N3 = L3 * a;
    N4 = L1 * b;
    N5 = L2 * b;
    N6 = L3 * b;

    N = [N1; N2; N3; N4; N5; N6];

    % derivatives wrt (xi, eta, mu)
    % dL1/dxi = -1, dL1/deta = -1
    % dL2/dxi =  1, dL2/deta =  0
    % dL3/dxi =  0, dL3/deta =  1
    %
    % da/dmu = -1/2, db/dmu = 1/2

    % Node 1: N1 = L1 * a
    dN1_dxi  = -a;
    dN1_deta = -a;
    dN1_dmu  = L1 * (-1/2);

    % Node 2: N2 = L2 * a
    dN2_dxi  =  a;
    dN2_deta =  0;
    dN2_dmu  = L2 * (-1/2);

    % Node 3: N3 = L3 * a
    dN3_dxi  =  0;
    dN3_deta =  a;
    dN3_dmu  = L3 * (-1/2);

    % Node 4: N4 = L1 * b
    dN4_dxi  = -b;
    dN4_deta = -b;
    dN4_dmu  = L1 * (1/2);

    % Node 5: N5 = L2 * b
    dN5_dxi  =  b;
    dN5_deta =  0;
    dN5_dmu  = L2 * (1/2);

    % Node 6: N6 = L3 * b
    dN6_dxi  =  0;
    dN6_deta =  b;
    dN6_dmu  = L3 * (1/2);

    % assemble B: each column is [dNi/dxi; dNi/deta; dNi/dmu], 6 by 3
    B = [ dN1_dxi  dN2_dxi  dN3_dxi  dN4_dxi  dN5_dxi  dN6_dxi;
          dN1_deta dN2_deta dN3_deta dN4_deta dN5_deta dN6_deta;
          dN1_dmu  dN2_dmu  dN3_dmu  dN4_dmu  dN5_dmu  dN6_dmu]';
end
