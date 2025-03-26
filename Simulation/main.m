%% Presim
clear
close all
clc



%% 基本参数

global e_previous sigma_previous j level W_stack W_optimal r_sum W_stack_unsafe W_optimal_unsafe
a1 = -1.3; A1 = 0.5; %Barrier Boundaries for x1, x1-coordinate values will be remained within (a1,A1) 
a2 = -3.1; A2 = 0.5; %Barrier Boundaries for x2: x2-coordinate values will be remained within (a2,A2) 

P.n = 2; %number of states
P.m = 1; %number of dimension
P.l = 4; %number of unknown parameters
P.W = 3; % number of the Actor Weight vector elements as wellas Critic Weight vector elements.
P.G =3;  % number of the Gamma vector elements. 

K = 5;

W_stack = zeros(2*(K+1),3);
W_optimal = zeros(2,3);

P.x0 = [-1;-3]; % Initial x-coordinate values, it has to be within (a,A).
% P.x0 = [-0.7;-2]; % Initial x-coordinate values, it has to be within (a,A).



%% 坐标变换

%transforming initial x-coordinates to initial s-coordinates. 
s10 = (log ((A1/a1)*((a1-P.x0(1))/(A1-P.x0(1)))));
s20 = (log ((A2/a2)*((a2-P.x0(2))/(A2-P.x0(2)))));

P.s0 = [s10;s20]; % Initial s-coordinate values

%% human-machine both safe

% Bounded Level-k Rationality 分层理性

options = odeset('OutputFcn',@odeplot);
a0 = 1.5;  % 3以上

for level = 0:K-1
    for j=0:1
        e_previous = zeros(2,10);
        sigma_previous = zeros(3,10);
        W0 = a0*ones(4*P.W,1);
        P.z0 = [P.s0;W0;P.x0]';
%         [t, z] = ode23('human_machine_safe',[0 30],P.z0,options);
        [t, z] = ode23('human_machine_safe',[0 30],P.z0);
        W_stack((level+1)*2+1+j,:) = z(end,3+6*j:5+6*j);
        figure
        plot(z(:,15),z(:,16),'-s','Linewidth',2)
        rectangle('position',[-1.3 -3.1 1.8 3.6])
        axis([-1.8 1 -3.6 1])
        hold off
    end
end

% Optimality Learning 理性层级 = ∞

j = 2;
e_previous = zeros(2,10);
sigma_previous = zeros(3,10);
P.z0 = [P.s0;a0*ones(4*P.W,1);P.x0]';
% [t, z] = ode23('closedLoopDynamics',[0 30],P.z0,options);
[t, z] = ode23('human_machine_safe',[0 30],P.z0);
figure
plot(z(:,15),z(:,16),'-s','Linewidth',2)
rectangle('position',[-1.3 -3.1 1.8 3.6])
axis([-1.8 1 -3.6 1])
hold off

W_optimal(1,:) = z(end,3:5);
W_optimal(2,:) = z(end,9:11);
W_stack((K+1)*2+1,:) = W_optimal(1,:);
W_stack((K+1)*2+2,:) = W_optimal(2,:);



%% human-machine both unsafe
% Bounded Level-k Rationality 分层理性
W_stack_unsafe = zeros(2*(K+1),3);
W_optimal_unsafe = zeros(2,3);
options = odeset('OutputFcn',@odeplot);
a0 = 1.5;  % 3以上

for level = 0:K-1
    for j=0:1
        e_previous = zeros(2,10);
        sigma_previous = zeros(3,10);
        W0 = a0*ones(4*P.W,1);
        P.z0 = [P.x0;W0]';
        [t, z] = ode23('human_machine_unsafe',[0 30],P.z0);
        W_stack_unsafe((level+1)*2+1+j,:) = z(end,3+6*j:5+6*j);
        figure
        plot(z(:,1),z(:,2),'-s','Linewidth',2)
        rectangle('position',[-1.3 -3.1 1.8 3.6])
        axis([-1.8 1 -3.6 1])
        hold off
    end
end

% Optimality Learning 理性层级 = ∞

j = 2;
e_previous = zeros(2,10);
sigma_previous = zeros(3,10);
P.z0 = [P.x0;a0*ones(4*P.W,1)]';
% [t, z] = ode23('closedLoopDynamics',[0 30],P.z0,options);
[t_machine_optimal, z_machine_optimal] = ode23('human_machine_unsafe',[0 30],P.z0);
figure
plot(z_machine_optimal(:,1),z_machine_optimal(:,2),'-s','Linewidth',2)
rectangle('position',[-1.3 -3.1 1.8 3.6])
axis([-1.8 1 -3.6 1])
hold off



%% 存储数据

W_optimal_unsafe(1,:) = z(end,3:5);
W_optimal_unsafe(2,:) = z(end,9:11);
W_stack_unsafe((K+1)*2+1,:) = W_optimal_unsafe(1,:);
W_stack_unsafe((K+1)*2+2,:) = W_optimal_unsafe(2,:);


%% bar，画模拟人行为概率分布图
% j = 4;
T = 0.1;
options = odeset('OutputFcn',@odeplot);
% r_sum = [];%
r_sum = zeros(1,7);
% lambda = [];
% lambda_sum = 0;

for j = 0:6
    x0=[P.s0;P.x0;0];
    for i = 1:300
        [t,x]= ode23('closedLoopDynamics_estimation_bar',[(i-1)*T i*T],x0);
        x0=[x(end,1:2)';x(end,3:4)';0];
        r_sum(j+1) = r_sum(j+1) + x(end,5);
    end
end


%% human-modeling，概率+仲裁函数
% j = 4;
T = 0.1;
options = odeset('OutputFcn',@odeplot);
% lambda = [];
% lambda_sum = 0;
% for j = 0:6
%     x0=[P.s0;P.x0];
%     for i = 1:300
%         [t,x]= ode23('closedLoopDynamics_human',[(i-1)*T i*T],x0);
%         x0=[x(end,1:2)';x(end,3:4)'];
%     end
% end

P.z0 = [P.s0;2*a0*ones(4*P.W,1);P.x0]';
[t_human_machine_online, z_human_machine_online] = ode23('closedLoopDynamics_human_machine_online',[0 30],P.z0);
figure
plot(z_human_machine_online(:,15),z_human_machine_online(:,16),'-s','Linewidth',2)
rectangle('position',[-1.3 -3.1 1.8 3.6])
axis([-1.8 1 -3.6 1])
hold off



P.z0 = [P.s0;P.x0]';
[t_only_human, z_only_human] = ode23('closedLoopDynamics_only_human',[0 30],P.z0);

P.z0 = [P.s0;P.x0]';
[t_only_human1, z_only_human1] = ode23('closedLoopDynamics_only_human1',[0 30],P.z0);

P.z0 = [P.s0;P.x0]';
[t_only_machine, z_only_machine] = ode23('closedLoopDynamics_only_machine',[0 30],P.z0);


figure
plot(t_human_machine_online,z_human_machine_online(:,15),LineWidth=2)
hold on
plot(t_only_human, z_only_human(:,3),LineWidth=2)
hold on
plot(t_only_human1, z_only_human1(:,3),LineWidth=2)
hold on
plot(t_machine_optimal, z_machine_optimal(:,1),LineWidth=2)
hold off
legend('human-machine','human','human1','machine')
% %% 子图像
% H = axes('Position',[0.18,0.62,0.28,0.25]); % 生成子图
% plot(t_only_human, z_only_human(:,3));                % 绘制局部曲线图
% xlim([min(t_only_human),max(t_only_human)]);    % 设置坐标轴范围                                                               
% set(H, 'XTick',[], 'YTick', []);



figure
plot(t_human_machine_online,z_human_machine_online(:,16),LineWidth=2)
hold on
plot(t_only_human, z_only_human(:,4),LineWidth=2)
hold on
plot(t_only_human1, z_only_human1(:,4),LineWidth=2)
hold on
plot(t_machine_optimal, z_machine_optimal(:,2),LineWidth=2)
hold off
legend('human-machine','human','human1','machine')
% 
% 
% figure
% plot(t_human_machine,z_only_machine(:,3)-z_human_machine(:,3),LineWidth=2)
% hold on
% plot(t_only_human, z_only_machine(:,3)-z_only_human(:,3),LineWidth=2)
% hold on
% plot(t_only_human1, z_only_machine(:,3)-z_only_human1(:,3),LineWidth=2)
% hold on
% plot(t_only_machine, z_only_machine(:,3)-z_only_machine(:,3),LineWidth=2)
% hold off
% legend('human-machine','human','human1','machine')
% 
% 




%% Online Learning of Bounded Level-k Rationality
% 
% options = odeset('OutputFcn',@odeplot);
% a0 = 1.5;  % 3以上
% 
% for level = 0:K-1
%     for j=0:1
%         e_previous = zeros(2,10);
%         sigma_previous = zeros(3,10);
%         W0 = a0*ones(4*P.W,1);
%         P.z0 = [P.s0;W0;P.x0]';
% %         [t, z] = ode23('human_machine_safe',[0 30],P.z0,options);
%         [t, z] = ode23('human_machine_safe',[0 30],P.z0);
% %         W_stack((level+1)*2+1+j,:) = z(end,3+6*j:5+6*j);
%         figure
%         plot(z(:,15),z(:,16),'-s','Linewidth',2)
%         rectangle('position',[-1.3 -3.1 1.8 3.6])
%         axis([-1.8 1 -3.6 1])
%         hold off
%     end
% end


