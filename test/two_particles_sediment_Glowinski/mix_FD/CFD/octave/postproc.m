% postproc.m: loads the data necessary and saves the graphs
% (compared to the reults by glowinski)

clear;
clc;
close all;

% read data from simulation
particle_pos_1 = load('../../DEM/post/position_particle_1.txt');
particle_pos_2 = load('../../DEM/post/position_particle_2.txt');

linienstaerke = 1;
MarkerGroesse = 4;

figure(1)
offset = particle_pos_1(1, 4);
h = plot(particle_pos_1(:, 1), particle_pos_1(:, 4), '*',
         particle_pos_2(:, 1), particle_pos_2(:, 4), '+');
set(h,'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(gca, 'FontSize', 12)
axis([0.0 0.7 0.0 0.04])
xlabel('time (s)')
ylabel('position (m)')
title('Comparison of the y-position of two particles', 'FontSize', 12)
legend('following particle', 'leading particle')
set(gca, 'FontSize', 12)
print('pos_y_two_particles.png')

clear;

% read data from simulation
particle_vel_1 = load('../../DEM/post/velocity_particle_1.txt');
particle_vel_2 = load('../../DEM/post/velocity_particle_2.txt');

linienstaerke = 1;
MarkerGroesse = 4;

figure(2)
h = plot(particle_vel_1(:, 1), particle_vel_1(:, 4), '-*',
         particle_vel_2(:, 1), particle_vel_2(:, 4), '-+');
set(h,'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(gca, 'FontSize', 12)
axis([0.0 0.7 -0.075 0.0])
xlabel('time (s)')
ylabel('z-veloctiy (m/s)')
title('Comparison of the settling velocity of two particles', 'FontSize', 12)
legend('following particle', 'leading particle')
set(gca, 'FontSize', 12)
print('vel_y_two_particles.png')

