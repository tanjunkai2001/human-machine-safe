%% Plot

% Phase portrait of the s-coordinates
figure(1)
plot(out.simout(:,1),out.simout(:,2),'-s','Linewidth',2)
xlabel('$s_{1}$(t)','interpreter','latex')
ylabel('$s_{2}$(t)','interpreter','latex')
set(gca,'FontSize',10)
l=legend('Transformed State Trajectory','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
% print -dpdf  States_RL_BF_FCL.pdf %saving plot as a pdf 
% save('States_RL_BF_FCL.mat') %saving data in the MAT format

% State trajectories of the s-coordinates
figure(2)
p=plot(out.tout,out.simout(:,1),out.tout,out.simout(:,2),'Linewidth',2);
mrk1={'s','v','o','*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time, t(s)','interpreter','latex')
ylabel('Transformed States, s (t)','interpreter','latex')
set(gca,'FontSize',10)
l=legend('$s_{1}$(t)','$s_{2}$(t)','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
% print -dpdf Time_Vs_TStates_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Time_Vs_TStates_RL_BF_FCL.mat') %saving data in the MAT format

% State trajectories of the x1-coordinates
figure(3)
p=plot(out.tout,out.simout(:,15),'-s','LineWidth',1);
xlabel('Time (s)','FontSize',16)
ylabel('$x_{1}$(t)','interpreter','latex')
set(gca,'FontSize',16)
hold on
p1=plot(out.tout,a1*ones(size(out.tout)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(out.tout,A1*ones(size(out.tout)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
l=legend('$x_{1}$ Trajectory','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
% print -dpdf  Time_Vs_Barrier_x1_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Time_Vs_Barrier_x1_RL_BF_FCL.mat') %saving data in the MAT format

% State trajectories of the x2-coordinates
figure(4)
p=plot(out.tout,out.simout(:,16),'-s','LineWidth',1);
xlabel('Time (s)','FontSize',16)
ylabel('$x_{2}$(t)','interpreter','latex')
set(gca,'FontSize',16)
hold on
p1=plot(out.tout,a1*ones(size(out.tout)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(out.tout,A2*ones(size(out.tout)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
l=legend('$x_{1}$ Trajectory','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
% print -dpdf  Time_Vs_Barrier_x2_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Time_Vs_Barrier_x2_RL_BF_FCL.mat') %saving data in the MAT format


% Trajectory of the Actor and Critic weights 
figure(5)
p=plot(out.tout,out.simout(:,3),out.tout,out.simout(:,4),out.tout,out.simout(:,5),out.tout,out.simout(:,6),out.tout,out.simout(:,7),out.tout,out.simout(:,8),'Linewidth',2);
mrk1={'s','v','o','*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time, t(s)','interpreter','latex')
ylabel('Weight Estimations for States, $\hat{W}$','interpreter','latex')
set(gca,'FontSize',10)
l=legend('$\hat{W_{a}}_{1}$','$\hat{W_{a}}_{2}$','$\hat{W_{a}}_{3}$','$\hat{W_{c}}_{1}$','$\hat{W_{c}}_{2}$','$\hat{W_{c}}_{3}$','interpreter','latex');
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
% print -dpdf  Actor_Weights_Estimations_for_States_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Actor_Weights_Estimations_for_States_RL_BF_FCL.mat') %saving data in the MAT format
 
% Phase portrait of the x-coordinates
figure(6)
plot(out.simout(:,15),out.simout(:,16),'-s','Linewidth',2)
rectangle('position',[-1.3 -3.1 1.8 3.6])
axis([-1.8 1 -3.6 1])
xlabel('$x_{1}$(t)','interpreter','latex')
ylabel('$x_{2}$(t)','interpreter','latex')
set(gca,'FontSize',10)
l=legend('Original State Trajectory','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
% print -dpdf  Main_States_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Main_States_RL_BF_FCL.mat') %saving data in the MAT format



% State trajectories of the x-coordinates
figure(7)
p=plot(out.tout,out.simout(:,15),out.tout,out.simout(:,16),'-s','Linewidth',2);
mrk1={'s','v','o','*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('Time, t(s)','interpreter','latex')
ylabel('Original States, s(t)','interpreter','latex')
set(gca,'FontSize',10)
l=legend('$x_{1}$(t)','$x_{2}$(t)','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
% print -dpdf  Time_Vs_OStates_RL_BF_FCL.pdf %saving plot as a pdf 
% save('Time_Vs_OStates_RL_BF_FCL.mat') %saving data in the MAT format
