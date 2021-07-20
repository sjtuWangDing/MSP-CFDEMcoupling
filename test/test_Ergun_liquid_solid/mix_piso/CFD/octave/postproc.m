close all;
clear;
clc;

%====================================%
% simulation data
%====================================%
path = '../postProcessing/volFlow_outlet/0/surfaceFieldValue_0.0005_B.dat';
data = load(path);
[nrow, ncol] = size(data)
% dp_sim = data(:, ncol - 1) - data(:, ncol) % 计算压力差
% t_sim = data(:, 1);

%====================================%
% properties data
%====================================%
rhoG = 1000            % density in kg/m3
dp = 0.0005            % particle diameter
phip = 1               % sphericity
Ustart = 0.0
Uend = 0.006
deltaU = (Uend - Ustart) / nrow
U = Ustart + deltaU:deltaU:Uend;
L = 0.0196             % length of bed
epsilon = 1.0 - (30000 * dp^3 * pi / 6) / (0.015^2 * pi * L / 4)  % voidfraction
nuG = 1.0*10^-6        % kinematic viscosity in m2/s
muG = nuG*rhoG         % dynamic viscosity in Pa*s
solidGravity = 2.65    % gravity of particle
outputArea = 1.7558150180*10^-4 % output face area

%===================
% Ergun Equation
%===================
dpErgun = (150 * ((1 - epsilon)^2 / epsilon^3) * ((muG .* U) / (phip * dp)^2)
            + 1.75 * ((1 - epsilon) / epsilon^3) * ((rhoG .* U.^2) / (phip*dp))
        );

%====================================%
% critical hydraulic gradient
%====================================%
mean_dp = 0.0005; % mean harmonic diameter
Zhou_beta = 3.5
Zhou_icr = (2.0 / 3.0) * (solidGravity - 1) * (dp^2) / (dp^2 + Zhou_beta * mean_dp^2 * epsilon^2 / (15 * (1 - epsilon)^2))
Terzaghi_icr = (1 - epsilon) * (solidGravity - 1)
Zhou_grad_p = Zhou_icr * (rhoG * 9.81)
Terzaghi_grad_p = Terzaghi_icr * (rhoG * 9.81)

%====================================%
% simulation data
%====================================%
sim_data = [
    [0.0, 0.0],
    [0.1, 7.64*10^-8],
    [0.2, 1.522*10^-7],
    [0.3, 2.273*10^-7],
    [0.4, 3.02*10^-7],
    [0.5, 3.762*10^-7],
    [0.6, 4.505*10^-7],
    [0.65, 4.873*10^-7],
    [0.7, 5.247*10^-7],
    [0.75, 5.618*10^-7],
    [0.8, 5.992*10^-7],
    [0.85, 6.366*10^-7],
    [0.9, 6.733*10^-7],
    [0.93, 7.015*10^-7],
    [0.939, 1.6*10^-6],
    [0.94, 2*10^-6],
]
Re = dp * 0.01 * rhoG / muG

%====================================%
% plot data
%====================================%

figure(1)
h0 = plot(sim_data(:, 1) * (rhoG * 9.81), sim_data(:, 2) / outputArea, '-s')
hold on
h1 = plot(dpErgun, U, '-bk')
hold on
axis([0, 10000, 0, 0.01])
ylim=get(gca,'Ylim');
h2 = plot([Zhou_grad_p, Zhou_grad_p], ylim, '--g');
h3 = plot([Terzaghi_grad_p, Terzaghi_grad_p], ylim, '--r');

linienstaerke = 1.5;
MarkerGroesse = 8;
set(h0, 'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(h1, 'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(h2, 'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(h3, 'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);

legend('Un-resolved (Numerical)', 'Eq.(xxx) (Ergun)', 'Eq.(xxx) (Zhou)', 'Eq.(xxx) (Terzaghi)', 'Location','North')
xlabel('Pressure gradient (Pa/m)')
ylabel('Output average velocity (m/s)')

print('superficial vel with hydraulic gradient.png')
