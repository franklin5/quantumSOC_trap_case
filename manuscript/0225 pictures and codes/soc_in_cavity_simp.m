
%Plot ep vs Omg
clear
clc

load('OmgcVsep_05to5.mat');
x=Expression1(1,:);
y1=Expression1(2,:);
y5=Expression1(6,:);
plot(x,y1,'k--','linewidth',2)
hold on
plot(x,y5,'r','linewidth',2)
set(gca,'fontsize',16)
xlabel('\epsilonp')
ylabel('\Omega')
axis([0 2.5 0 0.4])

