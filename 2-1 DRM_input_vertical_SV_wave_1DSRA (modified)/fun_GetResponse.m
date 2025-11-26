function [U,V] = fun_GetResponse(time,coords,type,angle)

global dt Disp Velo Acce
%roh = 1922.0;
Vs = 240;
poi = 1/3;
Vp = Vs*sqrt(2*(1-poi)/(1-2*poi));
theta_s = angle/180.0*pi;
theta_p = asin(Vp/Vs*sin(theta_s));
      
x0 = 9.5+10.5*tan(theta_s);
y0 = 0;
x_rela = coords(1) - x0;
y_rela = y0 - coords(2);
      
t_ini = 10.5/Vs/cos(theta_s);
%t_ini = 0;      
k = (Vp/Vs);
      
U_si = 1.0;
U_sr = (sin(2*theta_s)*sin(2*theta_p)-k^2*cos(2*theta_s)^2)/...
       (sin(2*theta_s)*sin(2*theta_p)+k^2*cos(2*theta_s)^2)*U_si;
U_pr = -(-2*k^2*sin(2*theta_s)*cos(2*theta_s))/...
       (sin(2*theta_s)*sin(2*theta_p)+k^2*cos(2*theta_s)^2)*U_si*(Vs/Vp);
  
t1 = -x_rela/Vs*sin(theta_s)+y_rela/Vs*cos(theta_s)+time-t_ini;
t2 = -x_rela/Vs*sin(theta_s)-y_rela/Vs*cos(theta_s)+time-t_ini;
t3 = -x_rela/Vp*sin(theta_p)-y_rela/Vp*cos(theta_p)+time-t_ini;

if type == 1
    t2 = 0;
    t3 = 0;
end
      
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
        
%     U(1) = U(1) + U_si*cos(theta_s)*interp1(t,Disp,t1);
%     U(2) = U(2) + U_si*cos(theta_s)*interp1(t,Velo,t1);
%     U(3) = U(3) + U_si*cos(theta_s)*interp1(t,Acce,t1);
%     
%     V(1) = V(1) + U_si*sin(theta_s)*interp1(t,Disp,t1);
%     V(2) = V(2) + U_si*sin(theta_s)*interp1(t,Velo,t1);
%     V(3) = V(3) + U_si*sin(theta_s)*interp1(t,Acce,t1);
    
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
        
%     U(1) = U(1) - U_sr*cos(theta_s)*interp1(t,Disp,t2);
%     U(2) = U(2) - U_sr*cos(theta_s)*interp1(t,Velo,t2);
%     U(3) = U(3) - U_sr*cos(theta_s)*interp1(t,Acce,t2);
%     
%     V(1) = V(1) + U_sr*sin(theta_s)*interp1(t,Disp,t2);
%     V(2) = V(2) + U_sr*sin(theta_s)*interp1(t,Velo,t2);
%     V(3) = V(3) + U_sr*sin(theta_s)*interp1(t,Acce,t2);
              
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
              
%     U(1) = U(1) + U_pr*sin(theta_p)*interp1(t,Disp,t3);
%     U(2) = U(2) + U_pr*sin(theta_p)*interp1(t,Velo,t3);
%     U(3) = U(3) + U_pr*sin(theta_p)*interp1(t,Acce,t3);
%     
%     V(1) = V(1) + U_pr*cos(theta_p)*interp1(t,Disp,t3);
%     V(2) = V(2) + U_pr*cos(theta_p)*interp1(t,Velo,t3);
%     V(3) = V(3) + U_pr*cos(theta_p)*interp1(t,Acce,t3);
end

V(1) = -V(1);
V(2) = -V(2);
V(3) = -V(3);

end

