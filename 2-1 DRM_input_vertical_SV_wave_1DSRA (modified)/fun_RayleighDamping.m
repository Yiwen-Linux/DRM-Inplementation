function aC = fun_RayleighDamping
%--------------------------------------------------------------------------
% Originally written by Wenyang Zhang
% Email: zwyll@ucla.edu
% Data : 09/12/2018
%
% Modified by Yiwen
% Data : 11/25/2025
%
% This code computes Rayleigh Damping coefficients
%----------------------------------------
% YOU MUST MODIFY THIS FUNCTION.
%
% Modification instruction:
%
%   - define rayleigh damping ratio xi1, xi2
%   - define first and second order frequency of soil column f1, f2
%   - you may directly output aC = [0,0] for an undamped system
%----------------------------------------
%
%
% INPUT:
%       (NONE)
%
% OUTPUT:
%       aC : 2 by 1 vecotr, Rayleigh Damping Coefficients.
%
%
%----------------------------------------
% YOU MUST MODIFY THIS FUNCTION.
%--------------------------------------------------------------------------

xi1 = 0.1;
xi2 = 0.1;
f1  = 0.35;
f2  = 20;
AC  =1/4/pi*[1/f1 4*pi^2*f1;1/f2 4*pi^2*f2];

aC = AC\[xi1;xi2];

%aC = [0;0]; undamped system

end