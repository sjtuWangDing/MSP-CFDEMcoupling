% postproc.m: loads the data necessary and saves the graphs
% (compared to the reults by glowinski)

clear;
clc;
close all;

% read data from simulation
particle_vel_rec = load('./timeAverageVel.txt');
exp_particle_vel = load('./timeAverageVel_exp.txt');

linienstaerke = 1;
MarkerGroesse = 6;

figure(1)
h = plot(particle_vel_rec(:, 4), particle_vel_rec(:, 1) - 0.026, '-or',
         exp_particle_vel(:, 1), exp_particle_vel(:, 2) / 1000.0, '-sbk');

set(h, 'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(gca, 'FontSize', 12)
axis([-0.05 0.25 0.0 0.08])
xlabel('Z component of time average particle veloctiy (m/s)')
ylabel('Height (m)')
legend('Numerical  ', 'Experimental  ')
print('time_ave_vel.png')
