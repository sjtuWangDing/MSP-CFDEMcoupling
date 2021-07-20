close all;
clear;
clc;

%====================================%
% simulation data
%====================================%
path = '../postProcessing/volFlow_outlet/0/surfaceFieldValue_semi.dat';
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
L = 0.0194             % length of bed
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
% Terzaghi_grad_p = Terzaghi_icr
% Zhou_grad_p = Zhou_icr

%====================================%
% simulation data
%====================================%
% sim_data_piso = [
%     [0.0, 0.0],
%     [0.1, 7.64*10^-8],
%     [0.2, 1.522*10^-7],
%     [0.3, 2.273*10^-7],
%     [0.4, 3.02*10^-7],
%     [0.5, 3.762*10^-7],
%     [0.6, 4.505*10^-7],
%     [0.65, 4.873*10^-7],
%     [0.7, 5.247*10^-7],
%     [0.75, 5.618*10^-7],
%     [0.8, 5.992*10^-7],
%     [0.85, 6.366*10^-7],
%     [0.9, 6.733*10^-7],
%     [0.93, 7.015*10^-7],
%     [0.94, 2*10^-6],
% ]

% sim_data_semi = [
%     [0.0, 0.0],
%     [0.1, 7.54*10^-8],
%     [0.2, 1.513*10^-7],
%     [0.3, 2.28*10^-7],
%     [0.4, 3.04*10^-7],
%     [0.5, 3.80*10^-7],
%     [0.6, 4.57*10^-7],
%     [0.65, 4.94*10^-7],
%     [0.7, 5.31*10^-7],
%     [0.75, 5.66*10^-7],
%     [0.8, 5.98*10^-7],
%     [0.85, 5.76*10^-7],
% ]

sim_data_piso = [
    [0.0, 0.0],
    [0.1, 7.28*10^-8],
    [0.2, 1.452*10^-7],
    [0.3, 2.17*10^-7],
    [0.4, 2.884*10^-7],
    [0.5, 3.595*10^-7],
    [0.6, 4.304*10^-7],
    [0.65, 4.657*10^-7],
    [0.7, 5.014*10^-7],
    [0.75, 5.367*10^-7],
    [0.8, 5.723*10^-7],
    [0.85, 6.011*10^-7],
    [0.9, 6.355*10^-7],
    [0.92, 6.493*10^-7],
    [0.94, 6.630*10^-7],
    [0.96, 6.767*10^-7],
    [0.97, 2*10^-6],
]

sim_data_semi = [
    [0.0, 0.0],
    [0.1, 7.198*10^-8],
    [0.2, 1.428*10^-7],
    [0.3, 2.129*10^-7],
    [0.4, 2.837*10^-7],
    [0.5, 3.537*10^-7],
    [0.6, 4.213*10^-7],
    [0.65, 4.589*10^-7],
    [0.7, 4.910*10^-7],
    [0.75, 5.245*10^-7],
    [0.8, 5.567*10^-7],
    [0.85, 5.870*10^-7],
    [0.9, 6.10*10^-7],
    [0.92, 6.40*10^-7],
    [0.94, 2*10^-6],
]
Re = dp * 0.01 * rhoG / muG

%====================================%
% plot data
%====================================%

figure(1)
h0 = plot(sim_data_piso(:, 1), sim_data_piso(:, 2) / outputArea, '-sr',
          sim_data_semi(:, 1), sim_data_semi(:, 2) / outputArea, '-^b',
          dpErgun / (rhoG * 9.81), U, '-bk'
)
axis([0, 1.0, 0, 0.01])
ylim=get(gca,'Ylim');
hold on
h1 = plot([Zhou_icr, Zhou_icr], ylim, '--m',
          [Terzaghi_icr, Terzaghi_icr], ylim, '--g'
);

linienstaerke = 1.5;
MarkerGroesse = 8;
set(h0, 'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(h1, 'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);

legend('Un-resolved (Numerical)', 'Cs-PFCM (Numerical)', 'Ergun (Theoretical)', 'Zhou (Theoretical)', 'Terzaghi (Theoretical)', 'Location','North')
xlabel('Pressure gradient, I')
ylabel('Output superficial velocity, Usup(m/s)')

print('superficial vel with hydraulic gradient.png')
