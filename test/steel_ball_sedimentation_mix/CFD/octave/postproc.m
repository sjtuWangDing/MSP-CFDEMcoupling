% postproc.m: loads the data necessary and saves the graphs (compared to the reults by glowinski)
%
close all;
clear;
clc;

%====================================%
% simulation data
%====================================%
path = '../../DEM/post/falling_height_velocity_particle.txt';
data = load(path);
X_sim = data(:,2:4);
U_sim = data(:,5:7);
t_sim = data(:,1);
U_mag = data(:,8);
height = 0.275; % 起始高度(m)
fprintf('final velocity of sim = %f/%f/%f m/s\n',U_sim(length(U_sim(:,1)),1),U_sim(length(U_sim(:,1)),2),U_sim(length(U_sim(:,1)),3) )

%====================================%
% plot data
%====================================%
plot(height-X_sim(:,3), U_mag, 'rd-*')
legend("simulation - Mix") 
xlabel('falling height (m)')
ylabel('velocity magnitude(m/s)')

%print('cfdemSolverPiso_settlingTestMPI.eps','-deps2')
print -color "cfdemSolverPISO.png"
