% postproc.m: loads the data necessary and saves the graphs
% (compared to the reults by glowinski)

clear;
clc;
close all;

% read data from simulation
particle_vel_1 = load('../../DEM/post/velocity_particle_1.txt');
linienstaerke = 1;
MarkerGroesse = 4;

figure(1)
h = plot(particle_vel_1(:, 1), particle_vel_1(:, 4), '*', 0.2, -0.08, '+');
set(h, 'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(gca, 'FontSize', 12)
axis([0.0 1.4 -0.18 0.0])
xlabel('time (s)')
ylabel('z-veloctiy (m/s)')
title('Velocity of the settling particle', 'FontSize', 12)
legend('particle')
print('vel_z_particles.png')
