% postproc.m: loads the data necessary and saves the graphs
% (compared to the reults by glowinski)

clear;
clc;
close all;

% read data from simulation
particle_vel_1 = load('./volAverageVel.txt');
linienstaerke = 1.5;
MarkerGroesse = 5;

nrow = size(particle_vel_1, 1)

figure(1)
h = plot(particle_vel_1(1:2:nrow+1, 1) - 0.95, particle_vel_1(1:2:nrow+1, 4), '-or');

xlim=get(gca,'Xlim');
hold on
plot(xlim, [0.182, 0.182], '--r');
plot(xlim, [0.172, 0.172], '--bk');

set(h, 'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(gca, 'FontSize', 12)
axis([0.0 10 0.05 0.3])
xlabel('time (s)')
ylabel('Z component of volume average particle veloctiy (m/s)')
legend('Numerical', 'Mean value(numerical)', 'Mean value(experimental   )')
print('vel_z_particles.png')
