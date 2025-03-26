figure
p=plot(t,z(:,15),'Linewidth',3)
hold on
p1=plot(t,a1*ones(size(t)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',2)
hasbehavior(p1, 'legend', false);
p1=plot(t,A1*ones(size(t)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',2)
hasbehavior(p1, 'legend', false);
axis([0 15 -4 2])
xlabel('t(s)','interpreter','latex')
ylabel('$x_{1}$(t)','interpreter','latex')
set(gca,'FontSize',20)
l=legend('Trajectory of state $x_{1}$(t)','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
title('State $x_{1}$(t)','interpreter','latex','FontSize',20,'fontname','Times New Roman')

% print -dpdf lcno-x3.pdf
% 导出到pdf
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
filename = 'state1'; % 设定导出文件名
print(gcf,filename,'-dpdf','-r0')
close(gcf)


figure
p=plot(t,z(:,16),'Linewidth',3)
hold on
p1=plot(t,a2*ones(size(t)),'--','LineWidth',1);
set(p1,'Color',get(p(1),'Color'),'LineWidth',1)
hasbehavior(p1, 'legend', false);
p1=plot(t,A2*ones(size(t)),'--','LineWidth',2);
set(p1,'Color',get(p(1),'Color'),'LineWidth',2)
hasbehavior(p1, 'legend', false);
axis([0 15 -4 2])
xlabel('t(s)','interpreter','latex')
ylabel('$x_{2}$(t)','interpreter','latex')
set(gca,'FontSize',20)
l=legend('Trajectory of state $x_{2}$(t)','interpreter','latex')
set(l,'Interpreter','latex','Location','northeast');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
title('State $x_{2}$(t)','interpreter','latex','FontSize',20,'fontname','Times New Roman')

% print -dpdf lcno-x3.pdf
% 导出到pdf
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
filename = 'state2'; % 设定导出文件名
print(gcf,filename,'-dpdf','-r0')
close(gcf)