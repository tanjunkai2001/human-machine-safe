figure
plot(z(:,15),z(:,16),'-s','Linewidth',2)
% rectangle('position',[-1.3 -3.1 1.8 3.6])
axis([-1.8 1 -3.6 1])
xlabel('$x_{1}$(t)','interpreter','latex')
ylabel('$x_{2}$(t)','interpreter','latex')
set(gca,'FontSize',20)
% l=legend('Original State Trajectory','interpreter','latex')
% set(l,'Interpreter','latex','Location','northeast');
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [6 5]);
% set(gcf, 'PaperPosition', [0 0 6 5]);
% grid on 



hold on

z_unsafe = load("unsafe.mat");
% figure
plot(z_unsafe.z(:,1),z_unsafe.z(:,2),'-s','Linewidth',2)
rectangle('position',[-1.3 -3.1 1.8 3.6],'Linewidth',2)
axis([-1.8 1 -3.6 1])
xlabel('$x_{1}$(t)','interpreter','latex')
ylabel('$x_{2}$(t)','interpreter','latex')
set(gca,'FontSize',20)

hold on

plot(0,0,'*','LineWidth',3)

title('State Trajectory','FontSize',20,'fontname','Times New Roman');

l=legend('Transformed-system Trajectory','Original-system Trajectory','Origin Point','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on

% % print -dpdf lcno-x3.pdf
% % 导出到pdf
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% filename = 'state_trajectory'; % 设定导出文件名
% print(gcf,filename,'-dpdf','-r0')
% close(gcf)