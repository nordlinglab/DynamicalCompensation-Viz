%% Rain200519_Phase_portrait_COVID-19_test_model - Article_Ashyani2019_Dynamical_compensation
% System:
%     dy/dt = by(t)+d(t)+sz(t)(lr(t)-y(t))
%     dz/dt = -cz(t)(r(t)-y(t))
clear all, clc, close all
% Initial setting
tf = 100;
tps = 0; % Pulse start time
dt = 0.01; %Sampling time
tpi = 0; % Pulse interval
tr = tf-tps-tpi; % Time points after pulse input
%% Figiure 2A
sys_a = 0; sys_d = 0.01; sys_r = 11; sys_l = 0.7;
sys_s = 0.25; sys_b = 0.3; sys_c = 2;
sys_E1 = [-sys_d/sys_b ; 0]; % equibrium point 1 [y; z]
sys_E2 = [sys_r; -(sys_b+sys_d/sys_r)*(1/(sys_s*(sys_l-1)))];

arrow_range = 11;
interval = 1;
color = cbrewer2('qual', 'Set1', 5);

% Solve ode
my_opts = odeset('RelTol',1e-8,'AbsTol',1e-10);  
f = @(t,x)DC_2D_generalized_model_phase(t,x,tps,tpi,sys_r,sys_r,...
    sys_a,sys_b,sys_d,sys_s,sys_l,sys_c); %[y; z]
figure(1), hold on
plot_initial_y = [9:0.25:10];
plot_initial_z = 6;

for j = plot_initial_y
    for s = plot_initial_z
        [ts,ys] = ode45(f, [0 tf], [j; s], my_opts);
        plot(ys(:,1),ys(:,2),'LineWidth',1, 'Color', color(2,:))
        hold on
    end
end

% Draw one of the trajectories as a red color
plot_initial_y_special = 10;
plot_initial_z_special = 6;
[ts,ys] = ode45(f, [0 tf], [j; s], my_opts);
plot(ys(:,1),ys(:,2),'LineWidth',1, 'Color', color(4,:), linewidth=2)

xlabel('$y(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
ylabel('$z(t)$', 'Fontsize', 12,'Interpreter', 'Latex')
set(gca,'FontSize',20)
grid on

% Define the range for phase plot
y1 = sys_E2(1)-arrow_range:interval:sys_E2(1):interval:sys_E2(1)+arrow_range;
y2 = sys_E2(2)-arrow_range:interval:sys_E2(2):interval:sys_E2(2)+arrow_range;
% Phase portrait
[x,y] = meshgrid(y1,y2);
mu = zeros(size(x));
v = zeros(size(x));
for i = 1:numel(x)
    Yprime = f(0,[x(i); y(i)]);
    mu(i) = Yprime(1);
    v(i) = Yprime(2);
end
quiver(x,y,mu,v,'Color', color(2,:), 'LineWidth', 0.5,...
    'AutoScaleFactor', 2);

figure(1), hold on
plot(sys_E1(1), sys_E1(2), '.', 'MarkerSize', 20, 'Color', color(1,:))
plot(sys_E2(1), sys_E2(2), '.', 'MarkerSize', 20, 'Color', color(1,:))
xlim([-0.2, 15])
ylim([-0.2, 15])

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 11, 16], 'PaperUnits', 'Inches', 'PaperSize', [10, 10])
saveas(gcf,'Akram_2023_DC_phase_portrait_original.pdf')
%% Figiure 2B
sys_a = 0; sys_d = 0.01; sys_r = 11; sys_l = 0.7;
sys_s = 1.5; sys_b = 0.3; sys_c = 2;
sys_E1 = [-sys_d/sys_b ; 0]; % equibrium point 1 [y; z]
sys_E2 = [sys_r; -(sys_b+sys_d/sys_r)*(1/(sys_s*(sys_l-1)))];

arrow_range = 14;
interval = 1;
color = cbrewer2('qual', 'Set1', 5);

% Solve ode
my_opts = odeset('RelTol',1e-8,'AbsTol',1e-10);  
f = @(t,x)DC_2D_generalized_model_phase(t,x,tps,tpi,sys_r,sys_r,...
    sys_a,sys_b,sys_d,sys_s,sys_l,sys_c); %[y; z]
figure(2), hold on
plot_initial_y = [9:0.25:10];
plot_initial_z = 6;

for j = plot_initial_y
    for s = plot_initial_z
        [ts,ys] = ode45(f, [0 tf], [j; s], my_opts);
        plot(ys(:,1),ys(:,2),'LineWidth',1, 'Color', color(2,:))
        hold on
    end
end

% Draw one of the trajectories as a red color
plot_initial_y_special = 10;
plot_initial_z_special = 6;
[ts,ys] = ode45(f, [0 tf], [j; s], my_opts);
plot(ys(:,1),ys(:,2),'LineWidth',1, 'Color', color(4,:), linewidth=2)

xlabel('$y(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
ylabel('$z(t)$', 'Fontsize', 12,'Interpreter', 'Latex')
set(gca,'FontSize',20)
grid on


% Define the range for phase plot
y1 = sys_E2(1)-arrow_range:interval:sys_E2(1):interval:sys_E2(1)+arrow_range;
y2 = sys_E2(2)-arrow_range:interval:sys_E2(2):interval:sys_E2(2)+arrow_range;
% Phase portrait
[x,y] = meshgrid(y1,y2);
mu = zeros(size(x));
v = zeros(size(x));
for i = 1:numel(x)
    Yprime = f(0,[x(i); y(i)]);
    mu(i) = Yprime(1);
    v(i) = Yprime(2);
end
quiver(x,y,mu,v,'Color', color(2,:), 'LineWidth', 0.5,...
    'AutoScaleFactor', 2);

figure(2), hold on
plot(sys_E1(1), sys_E1(2), '.', 'MarkerSize', 20, 'Color', color(1,:))
plot(sys_E2(1), sys_E2(2), '.', 'MarkerSize', 20, 'Color', color(1,:))
xlim([-0.2, 15])
ylim([-0.2, 15])

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 11, 16], 'PaperUnits', 'Inches', 'PaperSize', [10, 10])
saveas(gcf,'Akram_2023_DC_phase_portrait_change_s.pdf')
%% Figiure 2C
sys_a = 0; sys_d = 0.01; sys_r = 11; sys_l = 0.7;
sys_s = 0.25; sys_b = 0.6; sys_c = 2;
sys_E1 = [-sys_d/sys_b ; 0]; % equibrium point 1 [y; z]
sys_E2 = [sys_r; -(sys_b+sys_d/sys_r)*(1/(sys_s*(sys_l-1)))];

arrow_range = 10;
interval = 1;
color = cbrewer2('qual', 'Set1', 5);

% Solve ode
my_opts = odeset('RelTol',1e-8,'AbsTol',1e-10);  
f = @(t,x)DC_2D_generalized_model_phase(t,x,tps,tpi,sys_r,sys_r,...
    sys_a,sys_b,sys_d,sys_s,sys_l,sys_c); %[y; z]
figure(3), hold on
plot_initial_y = [9:0.25:10];
plot_initial_z = 6;

for j = plot_initial_y
    for s = plot_initial_z
        [ts,ys] = ode45(f, [0 tf], [j; s], my_opts);
        plot(ys(:,1),ys(:,2),'LineWidth',1, 'Color', color(2,:))
        hold on
    end
end

% Draw one of the trajectories as a red color
plot_initial_y_special = 10;
plot_initial_z_special = 6;
[ts,ys] = ode45(f, [0 tf], [j; s], my_opts);
plot(ys(:,1),ys(:,2),'LineWidth',1, 'Color', color(4,:), linewidth=2)

xlabel('$y(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
ylabel('$z(t)$', 'Fontsize', 12,'Interpreter', 'Latex')
set(gca,'FontSize',20)
grid on

% Define the range for phase plot
y1 = sys_E2(1)-arrow_range:interval:sys_E2(1):interval:sys_E2(1)+arrow_range;
y2 = sys_E2(2)-arrow_range:interval:sys_E2(2):interval:sys_E2(2)+arrow_range;
% Phase portrait
[x,y] = meshgrid(y1,y2);
mu = zeros(size(x));
v = zeros(size(x));
for i = 1:numel(x)
    Yprime = f(0,[x(i); y(i)]);
    mu(i) = Yprime(1);
    v(i) = Yprime(2);
end
quiver(x,y,mu,v,'Color', color(2,:), 'LineWidth', 0.5,...
    'AutoScaleFactor', 2);

figure(3), hold on
plot(sys_E1(1), sys_E1(2), '.', 'MarkerSize', 20, 'Color', color(1,:))
plot(sys_E2(1), sys_E2(2), '.', 'MarkerSize', 20, 'Color', color(1,:))

xlim([-0.2, 15])
ylim([-0.2, 15])


set(gcf, 'Units', 'Inches', 'Position', [0, 0, 11, 16], 'PaperUnits', 'Inches', 'PaperSize', [10, 10])
saveas(gcf,'Akram_2023_DC_phase_portrait_change_b.pdf')
%% Figiure 2D
sys_a = 0; sys_d = 0.01; sys_r = 11; sys_l = 0.7;
sys_s = 0.25; sys_b = 0.3; sys_c = 4;
sys_E1 = [-sys_d/sys_b ; 0]; % equibrium point 1 [y; z]
sys_E2 = [sys_r; -(sys_b+sys_d/sys_r)*(1/(sys_s*(sys_l-1)))];

arrow_range = 11;
interval = 1;
color = cbrewer2('qual', 'Set1', 5);

% Solve ode
my_opts = odeset('RelTol',1e-8,'AbsTol',1e-10);  
f = @(t,x)DC_2D_generalized_model_phase(t,x,tps,tpi,sys_r,sys_r,...
    sys_a,sys_b,sys_d,sys_s,sys_l,sys_c); %[y; z]
figure(4), hold on
plot_initial_y = [9:0.25:10];
plot_initial_z = 6;

for j = plot_initial_y
    for s = plot_initial_z
        [ts,ys] = ode45(f, [0 tf], [j; s], my_opts);
        plot(ys(:,1),ys(:,2),'LineWidth',1, 'Color', color(2,:))
        hold on
    end
end

% Draw one of the trajectories as a red color
plot_initial_y_special = 10;
plot_initial_z_special = 6;
[ts,ys] = ode45(f, [0 tf], [j; s], my_opts);
plot(ys(:,1),ys(:,2),'LineWidth',1, 'Color', color(4,:), linewidth=2)

xlabel('$y(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
ylabel('$z(t)$', 'Fontsize', 12,'Interpreter', 'Latex')
set(gca,'FontSize',20)
grid on

% Define the range for phase plot
y1 = sys_E2(1)-arrow_range:interval:sys_E2(1):interval:sys_E2(1)+arrow_range;
y2 = sys_E2(2)-arrow_range:interval:sys_E2(2):interval:sys_E2(2)+arrow_range;
% Phase portrait
[x,y] = meshgrid(y1,y2);
mu = zeros(size(x));
v = zeros(size(x));
for i = 1:numel(x)
    Yprime = f(0,[x(i); y(i)]);
    mu(i) = Yprime(1);
    v(i) = Yprime(2);
end
quiver(x,y,mu,v,'Color', color(2,:), 'LineWidth', 0.5,...
    'AutoScaleFactor', 2);

figure(4), hold on
plot(sys_E1(1), sys_E1(2), '.', 'MarkerSize', 20, 'Color', color(1,:))
plot(sys_E2(1), sys_E2(2), '.', 'MarkerSize', 20, 'Color', color(1,:))

xlim([-0.2, 15])
ylim([-0.2, 15])


set(gcf, 'Units', 'Inches', 'Position', [0, 0, 11, 16], 'PaperUnits', 'Inches', 'PaperSize', [10, 10])
saveas(gcf,'Akram_2023_DC_phase_portrait_change_c.pdf')
