function [U,V,W] = fun_GetResponse_3D(time,coords,angle,Vs,poi,Node_0,Surface_z)
%--------------------------------------------------------------------------
% Originally written by Wenyang Zhang
% Email: zwyll@ucla.edu
% Data : 09/12/2018
%
% This code generates analytical solutions for the free-field motions under
% incliend incident SV wave
%
% Inputs
%
% Time: current time
% coords: coordinates of the input node with format of [x, y, z]
% Vs: shear wave velocity of the homogeneous soil layer
% poi: Poisson's ratio of the homogeneous soil layer
% Node_0: the coordinates of the first node that will be excited with 
%         format of [x, y, z]. In most cases, this node belongs to one of 
%         four corner nodes on the bottom surface of the DRM interface
% Surface_z: the elevation of the ground surface. Here we assume z-axis 
%            represents the vertical direction.
%
%--------------------------------------------------------------------------
%
% Outputs
%
% U,V,W are three 1 by 3 vector that consist of the displacement, velocity,
% and acceleration values in x, y and z axis, respectively.
% E.g., U(1) = u_x, U(2) = v_x, U(3) = a_x
%
%--------------------------------------------------------------------------

global dt Disp Velo Acce

Vp = Vs*sqrt(2*(1-poi)/(1-2*poi));

alpha = angle(1);
beta  = angle(2);

theta_s = alpha/180.0*pi;
theta_p = asin(Vp/Vs*sin(theta_s));

theta_x = beta/180.0*pi;
      
x0 = Node_0(1) + (Surface_z-Node_0(3)) * tan(theta_s);
y0 = Surface_z;

n = [sin(theta_x),-cos(theta_x),0];
n_XoY = [cos(theta_x),sin(theta_x),0];

lambda = (coords-Node_0)*n';
Node_InPlane = coords-lambda*n;

if norm(norm(Node_InPlane-Node_0)) == 0
    n_InPlane = [0 0];
else
    n_InPlane = (Node_InPlane-Node_0)/norm(Node_InPlane-Node_0);

    if n_InPlane*n_XoY'>=0
        n_InPlane = [sqrt(1-n_InPlane(3)^2) n_InPlane(3)];
    else
        n_InPlane = [-sqrt(1-n_InPlane(3)^2) n_InPlane(3)];
    end
end

coords_InPlane = Node_0([1 3])+norm(Node_InPlane-Node_0)*n_InPlane;

x_rela = coords_InPlane(1) - x0;
y_rela = y0 - coords_InPlane(2);
      
t_ini = (Surface_z-Node_0(3))/Vs/cos(theta_s);

k = (Vp/Vs);
      
U_si = 1.0;
U_sr = (sin(2*theta_s)*sin(2*theta_p)-k^2*cos(2*theta_s)^2)/...
       (sin(2*theta_s)*sin(2*theta_p)+k^2*cos(2*theta_s)^2)*U_si;
U_pr = -(-2*k^2*sin(2*theta_s)*cos(2*theta_s))/...
       (sin(2*theta_s)*sin(2*theta_p)+k^2*cos(2*theta_s)^2)*U_si*(Vs/Vp);
  
t1 = -x_rela/Vs*sin(theta_s)+y_rela/Vs*cos(theta_s)+time-t_ini;
t2 = -x_rela/Vs*sin(theta_s)-y_rela/Vs*cos(theta_s)+time-t_ini;
t3 = -x_rela/Vp*sin(theta_p)-y_rela/Vp*cos(theta_p)+time-t_ini;
      
U(1) = 0.0;
U(2) = 0.0;
U(3) = 0.0;
V(1) = 0.0;
V(2) = 0.0;          
V(3) = 0.0;

if(t1 >= 0) 
    
    num1 = floor(t1/dt)+1;
    num2 = num1 + 1;
    t_ratio = (t1-(num1-1)*dt)/dt;
    
    U(1) = U(1) + U_si*cos(theta_s)*((1-t_ratio)*Disp(num1)+t_ratio*Disp(num2));
    U(2) = U(2) + U_si*cos(theta_s)*((1-t_ratio)*Velo(num1)+t_ratio*Velo(num2));
    U(3) = U(3) + U_si*cos(theta_s)*((1-t_ratio)*Acce(num1)+t_ratio*Acce(num2));
    
    V(1) = V(1) + U_si*sin(theta_s)*((1-t_ratio)*Disp(num1)+t_ratio*Disp(num2));
    V(2) = V(2) + U_si*sin(theta_s)*((1-t_ratio)*Velo(num1)+t_ratio*Velo(num2));
    V(3) = V(3) + U_si*sin(theta_s)*((1-t_ratio)*Acce(num1)+t_ratio*Acce(num2));

end
          
if(t2 >= 0)
    num1 = floor(t2/dt)+1;
    num2 = num1 + 1;
    t_ratio = (t2-(num1-1)*dt)/dt;
    
    U(1) = U(1) - U_sr*cos(theta_s)*((1-t_ratio)*Disp(num1)+t_ratio*Disp(num2));
    U(2) = U(2) - U_sr*cos(theta_s)*((1-t_ratio)*Velo(num1)+t_ratio*Velo(num2));
    U(3) = U(3) - U_sr*cos(theta_s)*((1-t_ratio)*Acce(num1)+t_ratio*Acce(num2));
    
    V(1) = V(1) + U_sr*sin(theta_s)*((1-t_ratio)*Disp(num1)+t_ratio*Disp(num2));
    V(2) = V(2) + U_sr*sin(theta_s)*((1-t_ratio)*Velo(num1)+t_ratio*Velo(num2));
    V(3) = V(3) + U_sr*sin(theta_s)*((1-t_ratio)*Acce(num1)+t_ratio*Acce(num2));
           
end
          
if(t3 >= 0) 
    num1 = floor(t3/dt)+1;
    num2 = num1 + 1;
    t_ratio = (t3-(num1-1)*dt)/dt;
    
    U(1) = U(1) + U_pr*sin(theta_p)*((1-t_ratio)*Disp(num1)+t_ratio*Disp(num2));
    U(2) = U(2) + U_pr*sin(theta_p)*((1-t_ratio)*Velo(num1)+t_ratio*Velo(num2));
    U(3) = U(3) + U_pr*sin(theta_p)*((1-t_ratio)*Acce(num1)+t_ratio*Acce(num2));
    
    V(1) = V(1) + U_pr*cos(theta_p)*((1-t_ratio)*Disp(num1)+t_ratio*Disp(num2));
    V(2) = V(2) + U_pr*cos(theta_p)*((1-t_ratio)*Velo(num1)+t_ratio*Velo(num2));
    V(3) = V(3) + U_pr*cos(theta_p)*((1-t_ratio)*Acce(num1)+t_ratio*Acce(num2));    
              
end

V(1) = -V(1);
V(2) = -V(2);
V(3) = -V(3);

W = V;
V = U*sin(theta_x);
U = U*cos(theta_x);

end

