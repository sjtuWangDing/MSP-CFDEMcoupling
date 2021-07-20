close all;
clear;
clc;

%====================================%
% simulation data
%====================================%
path = '../postProcessing/probes/0/p';
data = load(path);
[nrow, ncol] = size(data)
dp_sim = data(:, ncol - 1) - data(:, ncol); % 计算压力差
t_sim = data(:, 1);

%====================================%
% properties data
%====================================%
g = 9.81               % gravity m/s2
rhoG = 10              % density in kg/m3
rhoP = 2000            % particle density in kg/m3
dp = 0.001             % particle diameter
phip = 1               % sphericity
epsilon = 0.451335     % voidfraction
Ustart = 0.002         % start vel
Uend = 0.02            % end vel
timeStepNum = nrow
deltaU = (Uend - Ustart) / timeStepNum
U = Ustart + deltaU:deltaU:Uend
L = 0.0156             % length of bed
nuG = 1.5*10^-4        % kinematic viscosity in m2/s
muG = nuG*rhoG         % dynamic viscosity in Pa*s

%===================
% Ergun Equation
%===================
dpErgun= L * (150 * ((1 - epsilon)^2 / epsilon^3) * ((muG .* U) / (phip * dp)^2)
            + 1.75 * ((1 - epsilon) / epsilon^3) * ((rhoG .* U.^2) / (phip * dp))
        )

%==================================
% min fluidization velocity in m/s
%==================================
Umf = dp^2 * (rhoP - rhoG) * g / (150 * muG) * (epsilon^3 * phip^2) / (1 - epsilon)
ReMF = Umf * dp * rhoG / muG
if(ReMF<20)
    fprintf('applying eqn1 for Umf.\n')
elseif(ReMF>20 && ReMF<1000)
    fprintf('applying eqn1 for Umf.\n')
elseif (ReMF>=1000)
    fprintf('applying eqn2 for Umf.\n')
    Umf = sqrt(dp*(rhoP-rhoG)*g/(1.75*rhoG)*epsilon^3*phip);
    ReMF = Umf*dp*rhoG/muG;
end

dpUmf = L * (150 * ((1 - epsilon)^2 / epsilon^3) * ((muG .* Umf) / (phip * dp)^2)
            + 1.75 * ((1 - epsilon) / epsilon^3) * ((rhoG .* Umf.^2) / (phip * dp))
        )

%====================================%
% plot data
%====================================%
figure(1)
plot(U, dpErgun, '-', U, dp_sim * rhoG, '-', [Umf, Uend], dpUmf * ones(1,2), '--r')
print('superficial vel with hydraulic gradient.png')
