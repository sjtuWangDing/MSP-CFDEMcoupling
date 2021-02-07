% postproc.m: loads the data necessary and saves the graphs (compared to the reults by glowinski)

clear;

% read data from literature
load pos.mat

% read data from simulation
par1=load('../../DEM/post/position_particle_1.txt');
par2=load('../../DEM/post/position_particle_2.txt');

linienstaerke=2;
MarkerGroesse=8;

figure(1)
offset=5-par1(1,4);
h= plot(par1(:,1),offset+par1(:,4),'*',par2(:,1),offset+par2(:,4),'+',coord(1:dataLen(1,1),1+(1-1)*2),coord(1:dataLen(1,1),2+(1-1)*2),'o',coord(1:dataLen(2,1),1+(2-1)*2),coord(1:dataLen(2,1),2+(2-1)*2),'s');
 set(h,'LineWidth',linienstaerke,'MarkerSize',MarkerGroesse);
set(gca,'FontSize',14)
axis([0 .25 0 6])
xlabel('time (s)')
ylabel('position (cm)')
title('Comparison of the y-position of two particles','FontSize',15)
legend('following particle','leading particle','following particle Glowinski.','leading particle Glowinski.')
set(gca,'FontSize',12)
print('y-position.png')

clear;

% read data from literature
load vel.mat

% read data from simulation
par1vel=load('../../DEM/post/velocity_particle_1.txt');
par2vel=load('../../DEM/post/velocity_particle_2.txt');

linienstaerke=2;
MarkerGroesse=8;

figure(2)
h= plot(par1vel(:,1),par1vel(:,4),'*',par2vel(:,1),par2vel(:,4),'+',coord(1:dataLen(1,1),1+(1-1)*2),coord(1:dataLen(1,1),2+(1-1)*2),'o',coord(1:dataLen(2,1),1+(2-1)*2),coord(1:dataLen(2,1),2+(2-1)*2),'s');
 set(h,'LineWidth',linienstaerke,'MarkerSize',MarkerGroesse);
set(gca,'FontSize',14)
axis([0 .25 -20 0])
xlabel('time (s)')
ylabel('velocity (cm/s)')
title('Comparison of the settling velocity of two particles','FontSize',15)
legend('following particle','leading particle','following particle Glowinski.','leading particle Glowinski.')
set(gca,'FontSize',12)
print('y-velocity.png')
