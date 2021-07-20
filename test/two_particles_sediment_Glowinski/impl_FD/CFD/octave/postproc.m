% postproc.m: loads the data necessary and saves the graphs
% (compared to the reults by glowinski)

clear;
clc;
close all;

% read data from simulation
particle_pos_1 = load('./pos_1_simulation.txt');
particle_pos_2 = load('./pos_2_simulation.txt');

Apte_pos_1 = load('./pos_1_Apte.txt');
Apte_pos_2 = load('./pos_2_Apte.txt');

nrow_1 = size(particle_pos_1, 1)
nrow_2 = size(particle_pos_2, 1)
interval_nrow = 16;

linienstaerke = 1;
MarkerGroesse = 4;

figure(1)
offset = particle_pos_1(1, 4);
h = plot(Apte_pos_1(:, 1), Apte_pos_1(:, 2), '--^r',
         Apte_pos_2(:, 1), Apte_pos_2(:, 2), '-^r',
         particle_pos_1(2 : interval_nrow : nrow_1, 1), particle_pos_1(2 : interval_nrow : nrow_1, 4), '--sb',
         particle_pos_2(2 : interval_nrow : nrow_2, 1), particle_pos_2(2 : interval_nrow : nrow_2, 4), '-sb');
set(h,'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(gca, 'FontSize', 12)
axis([0.0 0.7 0.0 0.04])
xlabel('time (s)')
ylabel('position (m)')
title('Comparison of the z-position of two particles', 'FontSize', 12)
legend('P2, FDM (Apte)', 'P1, FDM (Apte)', 'P2, Cs-PFCM (this paper)', 'P1, Cs-PFCM (this paper)')
set(gca, 'FontSize', 12)
print('pos_z_two_particles.png')

clear;

% read data from simulation
particle_vel_1 = load('./vel_1_simulation.txt');
particle_vel_2 = load('./vel_2_simulation.txt');
vel_Idelsohn_1 = load('./vel_Idelsohn_1.txt');
vel_Idelsohn_2 = load('./vel_Idelsohn_2.txt');
vel_Apte_1 = load('./vel_Apte_1.txt');
vel_Apte_2 = load('./vel_Apte_2.txt');

linienstaerke = 1;
MarkerGroesse = 4;

nrow_1 = size(particle_vel_1, 1)
nrow_2 = size(particle_vel_2, 1)
interval_nrow = 4

% particle_vel_1(2 : interval_nrow : nrow_1, 1), particle_vel_1(2 : interval_nrow : nrow_1, 4), '-ob',
% particle_vel_2(2 : interval_nrow : nrow_2, 1), particle_vel_2(2 : interval_nrow : nrow_2, 4), '--ob',

figure(2)
h = plot(vel_Idelsohn_1(:, 1), vel_Idelsohn_1(:, 2), '-sm',
         vel_Idelsohn_2(:, 1), vel_Idelsohn_2(:, 2), '--sm',
         vel_Apte_1(:, 1), vel_Apte_1(:, 2), '-^r',
         vel_Apte_2(:, 1), vel_Apte_2(:, 2), '--^r',
         particle_vel_1(2 : interval_nrow : nrow_1, 1), particle_vel_1(2 : interval_nrow : nrow_1, 4), '-ob',
         particle_vel_2(2 : interval_nrow : nrow_2, 1), particle_vel_2(2 : interval_nrow : nrow_2, 4), '--ob'
);
set(h,'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(gca, 'FontSize', 12)
axis([0.0 0.7 -0.09 0.0])
xlabel('time (s)')
ylabel('z-veloctiy (m/s)')
title('Comparison of the settling velocity of two particles', 'FontSize', 12)
legend('P2, FDM (Idelsohn)', 'P1, FDM (Idelsohn)',
       'P2, FDM (Apte)', 'P1, FDM (Apte)',
       'P2, Cs-PFCM (this paper)', 'P1, Cs-PFCM (this paper)',
       'Location', 'North')
set(gca, 'FontSize', 12)
print('vel_z_two_particles.png')

