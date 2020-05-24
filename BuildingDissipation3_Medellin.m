function [ ] = BuildingDissipation3_Medellin()
%=======================================
%+++++++++++++++++++++++++++++++++++++++
%Inputs
%+++++++++++++++++++++++++++++++++++++++
%-----------------
%Von Mises System
%-----------------
Kvm=5;
m=1e-3;
%Fvm=10*sin(x*10*2*pi)
%Force in VM truss
x=linspace(-0.15,0.15,100);
Fvm=[];
for I=1:length(x); Fvm=[Fvm; VonMises(x(I))]; end
figure
plot(x,Fvm);
xlabel('Displacement [m]'); ylabel('Load [N]'); title('L-D Von Mises Reference')
print('FvmDual4','-dpng','-r300')
%--------------------------
%Building
%--------------------------
Kb=6.3;
Mb=1;
%--------------------------
%Dashpot
%--------------------------
eta=0.02;
natf=sqrt(Kb/Mb);
C=eta*2*natf*Mb;
C=C/2;
%--------------------------
%Force
%--------------------------
Medellin=load('SueloCMedellin.prn'); %Acceleration in gs
dt=0.02;
time=[0:dt:(length(Medellin)-1)*dt];
Medellin=[time' Medellin*9.8];

%+++++++++++++++++++++++++++++++++++++++
%Solution
%+++++++++++++++++++++++++++++++++++++++
t_start =0;
%t_end=1/f;
t_end=time(end);
time_span=linspace(t_start,t_end,1e6);

initial_position_M=0;
initial_speed_M=0;
initial_position_m=0;
initial_speed_m=0;

x0=[initial_position_M initial_speed_M initial_position_m initial_speed_m];

[t,x]=ode113(@(t,x) rhs(t,x,C,Kb,Mb,m,Kvm,Medellin),time_span,x0);

x1=x(:,1);
x1dot=x(:,2);
x2=x(:,3);
x2dot=x(:,4);
%Force Von Mises
Fvm=[];
for I=1:length(t); Fvm=[Fvm; VonMises(x1(I)-x2(I))]; end

Data=load('DataRefMedellin.dat');

figure
plot(Data(:,1),Data(:,2),'--k','LineWidth',1); hold on;
plot(t,x1,'-r','LineWidth',1.2); xlabel('Time [s]'); ylabel('Displacement [m]'); legend('Ref','Snap')
print('TimeDispX1_Dual4','-dpng','-r300')

figure
plot(Data(:,2),Data(:,3),'--k'); hold on;
plot(x1,Kb*x1+C*x1dot+Kvm*x2+C*x2dot,'-r'); xlabel('Displacement [m]'); ylabel('Load [N]');legend('Ref','Snap')
title('Reaction Left Wall')
print('FWall_Dual4','-dpng','-r300')

% dlmwrite('Case1.dat',[t x1 Kb*x1+C*x1dot+Kvm*x2])

figure
plot(x1-x2,Fvm,'-b','LineWidth',2); xlabel('Displacement VM [m]'); ylabel('Load [N]');
title('Load in Von Mises Truss -Simulation-');
hold on
% x=linspace(-0.5,0.5,100);
% Fvm=[];
% for I=1:length(x); Fvm=[Fvm; VonMises(x(I))]; end
% plot(x,Fvm,'--r'); legend('Simulation','Input')
% print('FVM_Dual4','-dpng','-r300')

figure
plot(x1,Kvm*x2,'--b'); hold on;
plot(x1,Fvm,'--r')
% plot(x1,Kb*x1,'--r');
% plot(x1,C*x1dot,'-m')
xlabel('Displacement X_1 [m]')
ylabel('Load [N]');
title('Force in Kvnm');
legend('Kvm','VM','Location','NorthWest')
legend box off
%print('Components4','-dpng','-r300')

figure
plot(t,x1-x2);
xlabel('t[s]'), ylabel('Displacement x_1-x_2 [m]')
title('Displacement Von Mises');
print('-dpng','-r300','DVM4')

max(abs(Data(:,2)))
max(abs(x1))

end
%----------------------------------------------
function [F]=VonMises(x)

k=12.537;
if x<=-0.1
    F=k*x+k*0.1;    
elseif x>=0.1
    F=k*x-k*0.1;
elseif x>=0 && x<0.1
    F=0.2*sin(x*10*2*pi);
else
    F=0.2*sin(x*10*2*pi);
end
end
%----------------------------------------------
function q=rhs(t,x,C,Kb,Mb,m,Kvm,Medellin)

aq1 = interp1(Medellin(:,1),Medellin(:,2),t);
%aq1 = interp1(Medellin(:,1),Medellin(:,2),tq,'spline');
%fprintf('t= %g, a=%g\n',t, aq1);

x1=x(1);
x1dot=x(2);
x2=x(3);
x2dot=x(4);

Fvm = VonMises(x1-x2);

%Dashpot parallel to kvn
x1dotdot=(Mb*aq1-C*x1dot-Kb*x1-Fvm)/Mb;
x2dotdot=(-Kvm*x2+Fvm-C*(x2dot))/m;

q=[x1dot x1dotdot x2dot x2dotdot]';

end
        
