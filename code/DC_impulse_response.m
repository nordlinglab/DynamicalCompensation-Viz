%% Noisy input and disturbance
% The system is defined by the following differential equations:
%     dy/dt = by(t)+d(t)-sz(t)(lr(t)-y(t))
%     dz/dt = -cz(t)(r(t)-y(t))
clear all, close all, clc

% Define color for the plot
color = cbrewer2('qual', 'Set1', 5);

% Set the seed for random number generator
rng(1);

%% Create input and disturbance
% Create noisy input
r = 11; % Define the value for r
sys_d = 0.01; % Define the system disturbance

input_shift = 50;
r_step = [ones(1, input_shift)*r, ones(1, 50)*r, ones(1, 250)*(r+0.5*r), ones(1, 50)*(r+0.25*r)];
r_noise = [zeros(1, input_shift), randn(1, 50), randn(1, 50), zeros(1, 150), randn(1, 50), zeros(1, 50)];
r_input = r_step + r_noise; % Add noise to the input
d_step = [ones(1, input_shift)*sys_d, ones(1, 200)*sys_d, ones(1, 100)*10, ones(1, 50)*5];
d_noise = [zeros(1, input_shift), zeros(1, 150), randn(1, 150), zeros(1, 50)];
d_input = d_step + d_noise; % Add noise to the disturbance
d_input(d_input < 0) = 0; % Ensure the disturbance is non-negative
tf = length(d_input); % Define the final time
t = 0:0.1:tf; % Define the time vector

%% Plot input
% Plot the input and disturbance
figure(200)
plot_input(r_input, d_input, t) % Call the function to plot the input and disturbance
xlabel('Time', 'Fontsize', 12, 'Interpreter', 'Latex') % Label the x-axis

% Set the figure size and save it as a PDF
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', 'Inches', 'PaperSize', [9.5, 9.5])
saveas(gcf,'Akram_2023_DC_noisy_step_input.pdf')

%% Case 1: Changing s
% Define the system parameters
sys_a = 0;  sys_l = 0.7; sys_s = 0.25; sys_b = 0.3; sys_c = 2;
parm1 = [sys_a, sys_l, sys_s, sys_b, sys_c]; % Store the parameters in a vector
x0 = [r; -(sys_b+sys_d/r)*(1/(sys_s*(sys_l-1)))]; % Define the initial conditions
x0_shift = [10, 6]; % Define the shift for the initial conditions

% Plot the input
figure(1)
plot_input(r_input, d_input, t)

% Plot the response comparison for changing s
sys_a = 0; sys_l = 0.7; sys_s = 0.25*6; sys_b = 0.3; sys_c = 2;
parm2 = [sys_a, sys_l, sys_s, sys_b, sys_c]; % Store the new parameters in a vector
plot_response(parm1, parm2, x0_shift, tf, r_input, d_input, color) % Call the function to plot the response

% Set the figure size and save it as a PDF
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 12], 'PaperUnits', 'Inches', 'PaperSize', [10-1.5, 12-1.5])
saveas(gcf,'Akram_2023_DC_noisy_step_input_change_s.pdf')

%% Case 2: Changing b
% Define the system parameters
sys_a = 0;  sys_l = 0.7; sys_s = 0.25; sys_b = 0.3; sys_c = 2;
parm1 = [sys_a, sys_l, sys_s, sys_b, sys_c]; % Store the parameters in a vector
x0 = [r; -(sys_b+sys_d/r)*(1/(sys_s*(sys_l-1)))]; % Define the initial conditions

% Plot the input
figure(2)
plot_input(r_input, d_input, t)

% Plot the response comparison for changing b
sys_a = 0; sys_d = 0.01; sys_l = 0.7; sys_s = 0.25; sys_b = 0.3*2; sys_c = 2;
parm2 = [sys_a, sys_l, sys_s, sys_b, sys_c]; % Store the new parameters in a vector
plot_response(parm1, parm2, x0_shift, tf, r_input, d_input, color) % Call the function to plot the response

% Set the figure size and save it as a PDF
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 12], 'PaperUnits', 'Inches', 'PaperSize', [10-1.5, 12-1.5])
saveas(gcf,'Akram_2023_DC_noisy_step_input_change_b.pdf')

%% Case 3: Changing c
% Define the system parameters
sys_a = 0;  sys_l = 0.7; sys_s = 0.25; sys_b = 0.3; sys_c = 2;
parm1 = [sys_a, sys_l, sys_s, sys_b, sys_c]; % Store the parameters in a vector
x0 = [r; -(sys_b+sys_d/r)*(1/(sys_s*(sys_l-1)))]; % Define the initial conditions

% Plot the input
figure(3)
plot_input(r_input, d_input, t)

% Plot the response comparison for changing c
sys_a = 0; sys_d = 0.01; sys_l = 0.7; sys_s = 0.25; sys_b = 0.3; sys_c = 2*2;
parm2 = [sys_a, sys_l, sys_s, sys_b, sys_c]; % Store the new parameters in a vector
plot_response(parm1, parm2, x0_shift, tf, r_input, d_input, color) % Call the function to plot the response

% Set the figure size and save it as a PDF
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 12], 'PaperUnits', 'Inches', 'PaperSize', [10-1.5, 12-1.5])
saveas(gcf,'Akram_2023_DC_noisy_step_input_change_c.pdf')

%% Function to plot input with noise
function plot_input(r_input, d_input, t)
    % This function plots the input and disturbance

    % Plot the input
    subplot(4,1,1)
    u = interp1([0:length(r_input)-1], r_input, t); % Interpolate the input
    plot(t, u, '-k', 'LineWidth', 2), hold on % Plot the input
    ylabel('$r(t)$', 'Fontsize', 12, 'Interpreter', 'Latex') % Label the y-axis
    set(gca,'FontSize',12)
    grid on
    
    % Plot the disturbance
    subplot(4,1,2)
    d = interp1([0:length(d_input)-1], d_input, t); % Interpolate the disturbance
    plot(t, d, '-k', 'LineWidth', 2), hold on % Plot the disturbance
    ylabel('$d(t)$', 'Fontsize', 12, 'Interpreter', 'Latex') % Label the y-axis
    set(gca,'FontSize',12)
    grid on
end

%% Plot state reponses of y and z
function plot_response(parm1, parm2, x0_shift, tf, r_input, d_input, color)
    sys_a = parm1(1); sys_l = parm1(2); sys_s = parm1(3); sys_b = parm1(4); sys_c = parm1(5);
    sys_b_case1 = sys_b; sys_s_case1 = sys_s; sys_c_case1 = sys_c;
    % ODE solver
    my_opts = odeset('RelTol',1e-8,'AbsTol',1e-10);  
    f = @(t,x)DC_2D_generalized_model(t,x,r_input,...
        sys_a,sys_b,d_input,sys_s,sys_l,sys_c); %[y; z]
    [t,y] = ode45(f, [0:0.8:tf], x0_shift, my_opts);
    y_old = y;
    
    subplot(4,1,3)
    p1 = plot(t, y(:,1), '-', 'LineWidth', 2, 'Color', color(2,:)), hold on
    ylabel('$y(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
    set(gca,'FontSize',12)
    grid on
    
    subplot(4,1,4)
    plot(t, y(:,2), '-', 'LineWidth', 2, 'Color', color(2,:)), hold on
    xlabel('Time', 'Fontsize', 12, 'Interpreter', 'Latex')
    ylabel('$z(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
    set(gca,'FontSize',12)
    grid on
    hold on

    sys_a = parm2(1); sys_l = parm2(2); sys_s = parm2(3); sys_b = parm2(4); sys_c = parm2(5);
    sys_b_case2 = sys_b; sys_s_case2 = sys_s; sys_c_case2 = sys_c;
    % ODE solver
    my_opts = odeset('RelTol',1e-8,'AbsTol',1e-10);  
    f = @(t,x)DC_2D_generalized_model(t,x,r_input,...
        sys_a,sys_b,d_input,sys_s,sys_l,sys_c); %[y; z]
    [t,y] = ode45(f, [0:0.8:tf], x0_shift, my_opts);
    difference = y-y_old;
    
    subplot(4,1,3)
    p2 = plot(t, y(:,1), '--', 'LineWidth', 2, 'Color', color(1,:)), hold on
    p3 = plot(t, difference(:,1), '-.', 'LineWidth', 2, 'Color', color(3,:))
    ylabel('$y(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
    set(gca,'FontSize',12)
    grid on
    ylim([-7, 24])
    legend([p3], 'Difference')

    subplot(4,1,4)
    plot(t, y(:,2), '--', 'LineWidth', 2, 'Color', color(1,:)), hold on
    xlabel('Time', 'Fontsize', 12, 'Interpreter', 'Latex')
    ylabel('$z(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
    set(gca,'FontSize',12)
    grid on
    if sys_b_case1 ~= sys_b_case2
        legend(['b=', num2str(sys_b_case1)], ['b=', num2str(sys_b_case2)])
    elseif sys_s_case1 ~= sys_s_case2
        legend(['s=', num2str(sys_s_case1)], ['s=', num2str(sys_s_case2)])
    elseif sys_c_case1 ~= sys_c_case2
        legend(['c=', num2str(sys_c_case1)], ['c=', num2str(sys_c_case2)])
    else
        disp("Parameter setting error")
    end
    
    ylim([0, 20])
    
end