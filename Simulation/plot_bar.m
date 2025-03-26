% %% bar
% % j = 4;
% T = 0.1;
% options = odeset('OutputFcn',@odeplot);
% % r_sum = [];%
% r_sum = zeros(1,7);
% lambda = [];
% lambda_sum = 0;
% % for i = 1:300
% %     x0=[P.s0;P.x0;zeros(5,1)];
% % %     [t,x]= ode23('dynamic3_estimation',[(i-1)*T i*T],x0,options);
% %     [t,x]= ode23('closedLoopDynamics_estimation_bar',[(i-1)*T i*T],x0);
% %     r_sum = r_sum + x(5:9);
% %     sigma = softmax(-x(5:9)');
% %     lambda_sum =  sigma'*[1:5]'*i/(i+1);
% %     lambda = [lambda lambda_sum];
% % end
% for j = 0:6
%     x0=[P.s0;P.x0;0];
%     for i = 1:300
%         [t,x]= ode23('closedLoopDynamics_estimation_bar',[(i-1)*T i*T],x0);
%         x0=[x(end,1:2)';x(end,3:4)';0];
%         r_sum(j+1) = r_sum(j+1) + x(end,5);
% %         r_sum = [r_sum;x(5)];
% %         lambda_sum =  sigma'*[1:5]'*i/(i+1);
% %         lambda = [lambda lambda_sum];
%     end
% end
chosen = [r_sum(1) r_sum(3) r_sum(4) r_sum(5) r_sum(6)];
% sigma = softmax(-r_sum');
sigma = softmax(-chosen');
% chosen = [sigma(1);sigma(3);sigma(4);sigma(5);sigma(6)];

figure
bar(sigma)
xlabel('Bounded Rationality Level','fontsize',20,'fontname','Times New Roman')
ylabel('$\mathcal{P}(\mathcal{U}_h = \hat{u}_h^k)$','Interpreter','latex','fontsize',20,'fontname','Times New Roman')
set(gca,'FontSize',16)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on 
title('Probabilistic Human Behavior Modeling','fontsize',20,'fontname','Times New Roman')