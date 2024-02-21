% This function models a 2D dynamical system governed by the following
% equations:
%   dy/dt = a + b*y + d + s*z*(l*r - y)
%   dz/dt = -c*z*(r-y)
%
% Inputs:
%   t: Time
%   x: States [y; z] or [y, z]
%   r_input: Reference input of y and time when it changes
%   a, b, d, s, l, c: Parameters of the system
%
% Outputs:
%   dydt: Change in y over time
%   dzdt: Change in z over time


function [dydt dzdt] = DC_2D_generalized_model(t,x,r_input,a,b,d_input,s,l,c)
% Input setting
r = interp1([0:length(r_input)-1],r_input,t);
d = interp1([0:length(d_input)-1],d_input,t);
% Block for caluclating may gradients at once
if prod(size(x)) == 2 && nargout == 1 %Single point give and only one output requested
    dydt = gradientcalc; % Preforms calculation of gradient in the single point
    return
else
    % Assumes that the function is used to report the gradient at the time
    % point t for many points, e.g. for plotting of phase diagram
    if size(x,1) == 2*size(x,2) % Assumes that the points are stacked vertically
        y = x(1:size(x,2),:);r
        z = x(size(x,2)+1:end,:);
    elseif size(x,2) == 2*size(x,1) % Assumes that the points are stacked horisontally
        y = x(:,1:size(x,1));
        z = x(:,size(x,1)+1:end);
    elseif size(x,3) == 2 && size(x,2) == size(x,1)
        y = x(:,:,1);
        z = x(:,:,2);
    else
        error('Unsupported format')
    end
    whos x y z
    x = zeros(2,1);
    py = NaN.*ones(size(y));
    pz = NaN.*ones(size(z));
    for i = 1:size(y,1)
        for j = 1:size(y,2)
            x(1) = y(i,j);
            x(2) = z(i,j);
            dxdt = gradientcalc; % Preforms calculation of gradient in the single point
            py(i,j) = dxdt(1);
            pz(i,j) = dxdt(2);
        end
    end
    dydt = py;
    dzdt = pz;
    return
end
        
function dxdt = gradientcalc
% Preforms calculation of gradient in a single point using access to
% internal variables in parent function
    % Adaptive proportional-integral feedback
    dxdt = [a+b*x(1)+d+s*x(2)*(l*r-x(1));
            -c*x(2)*(r-x(1))];
end
end