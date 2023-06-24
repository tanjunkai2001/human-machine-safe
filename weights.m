
figure
p=plot(t,z(:,6),t,z(:,7),t,z(:,8),'Linewidth',2);
mrk1={'s','v','o'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('t(s)','interpreter','latex')
ylabel('$\hat{W}_h$','interpreter','latex')

title('Critic Weights of Modeled Human','FontSize',20,'fontname','Times New Roman')
set(gca,'FontSize',20)
l=legend('$\hat{W}_{h,1}$','$\hat{W}_{h,2}$','$\hat{W}_{h,3}$','interpreter','latex');
set(l,'Interpreter','latex','Location','best');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 

% print -dpdf lcno-x3.pdf
% 导出到pdf
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
filename = 'human'; % 设定导出文件名
print(gcf,filename,'-dpdf','-r0')
close(gcf)


figure
p=plot(t,z(:,12),t,z(:,13),t,z(:,14),'Linewidth',2);
mrk1={'*','^','x'};
mrk=(mrk1(1,1:size(p,1)))';
set(p,{'marker'},mrk,{'markerfacecolor'},get(p,'Color'),'markersize',5);
f_nummarkers(p,10,1);
for i=1:size(p,1)
    hasbehavior(p(i), 'legend', false);
end
xlabel('t(s)','interpreter','latex')
ylabel('$\hat{W}_m$','interpreter','latex')
title('Critic Weights of the machine','FontSize',20,'fontname','Times New Roman')

set(gca,'FontSize',20)
l=legend('$\hat{W}_{m,1}$','$\hat{W}_{m,2}$','$\hat{W}_{m,3}$','interpreter','latex');
set(l,'Interpreter','latex','Location','best');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 

% print -dpdf lcno-x3.pdf
% 导出到pdf
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
filename = 'machine'; % 设定导出文件名
print(gcf,filename,'-dpdf','-r0')
close(gcf)

