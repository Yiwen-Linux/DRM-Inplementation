function [Tri_pts_array, Tri_wgt_array] = fun_Order2TriangleGaussQuadDataOrder
% 3 point Gaussian Quadrature on triangle (0,0)-(1,0)-(0,1)

    % (xi, eta)
    Tri_pts_array = [ 1/6, 1/6;
                 2/3, 1/6;
                 1/6, 2/3 ];   

    % weight
    Tri_wgt_array = [1/6;
               1/6;
               1/6];
end