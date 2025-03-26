function [zDot] = closedLoopDynamics_human_machine_online(t,z)
global e_previous sigma_previous j level W_stack W_optimal r_sum W_stack_unsafe W_optimal_unsafe
% Control Penalty Matrix  
R1=0.1*10;
R2=R1/2;
% State Penalty Matrix  
Q=[10 0; 0 10]/10; 
% Tuning Gains
etac1 = 50;
etac2 = 50;
etaa1 = 30;
etaa2 = 30;
beta = 30;
v = 1;

single = 1;

M= 10;
%% 
a1 = -1.3; A1 = 0.5; %Barrier Boundaries for sii(1), sii(1)-coordinate values will be remained within (a1,A1) 
a2 = -3.1; A2 = 0.5; %Barrier Boundaries for sii(2), sii(2)-coordinate values will be remained within (a2,A2) 

s = z(1:2,1); % row 1-2 in z vector
x1 = z(15,1);
x2 = z(16,1);
WaH1 = z(3:5); % row 35-37 in z vector
WcH1 = z(6:8); % row 38-40 in z vector
WaH2 = z(9:11); % row 35-37 in z vector
WcH2 = z(12:14); % row 38-40 in z vector

%% Dynamics

%%transforming initial s-coordinates to initial x-coordinates. 
x1 = a1*A1*((exp(s(1)))-1)/((a1*exp(s(1)))-A1);
x2 = a2*A2*((exp(s(2)))-1)/((a2*exp(s(2)))-A2);


%% Original dynamics without the parameters
F1 = [x2; -x2 - 0.5*x1 + 0.25*x2*(cos(2*x1)+2)^2 + 0.25*x2*(sin(2*x1)+2)^2];
G1 = [0; (cos(2*x1)+2)];
G2 = [0; (sin(4*x1^2)+2)];

%% Original dynamics without the parameters
F1x = [x2; -x2 - 0.5*x1 + 0.25*x2*(cos(2*x1)+2)^2 + 0.25*x2*(sin(2*x1)+2)^2];
G1x = [0; (cos(2*x1)-2)];
G2x = [0; (sin(4*x1^2)-2)];
%Known Parameters of the dynamics

F1_1x = F1x(1,:);
F1_2x = F1x(2,:);

G1_x = (cos(2*x1)+2);
G2_x = (sin(4*x1^2)+2);

F1_x = [F1_1x; F1_2x];
G1_x = [0; G1_x];
G2_x = [0; G2_x];

%% Terms in weight update laws that are evaluated along the system trajectory x

% Basis vector: [x1^2, x1x2,x2^2]
phix = [x1^2 x1*x2 x2^2];
phi_px = [2*x1 0; x2 x1; 0 2*x2]; % Jacobian of basis vector
%phi_p_ = [2*x1 0; x2 x1; 0 2*x2]; % Jacobian of basis vector

%% 更新律
%controller of the transformed system

mu1c = -0.5.*(R1\G1_x')*phi_px'*WcH1;
mu2c = -0.5.*(R2\G2_x')*phi_px'*WcH2;

%% Cost Function in the s-coordinate
x = [z(15,1);z(16,1)];
r1 = x'*Q*x + mu1c'*R1*mu1c + mu2c'*R2*mu2c;
r2 = r1/2;

%% Definition
omega = phi_px*(F1_x+G1_x*mu1c+G2_x*mu2c);
rho = (1+v*(omega'*omega))^2;

delta1x = WcH1'*omega + r1;
delta2x = WcH2'*omega + r2;

WcH1_ = -etac1*omega*delta1x/rho;% + clWc1; % WcH
WcH2_ = -etac2*omega*delta2x/rho;% + clWc2; % WcH

WaH1_ = - beta*(WaH1 - WcH1);
WaH2_ = - beta*(WaH2 - WcH2);

%% BE Extrapolation - Terms in weight update laws that are evaluated along arbitrarily selected trajectories sii

for i=1:M
    WcH1_ = WcH1_ - etac1*(sigma_previous(:,i)./(sigma_previous(:,i)'*sigma_previous(:,i)+1)^2)*e_previous(1,i)';
    WcH2_ = WcH2_ - etac2*(sigma_previous(:,i)./(sigma_previous(:,i)'*sigma_previous(:,i)+1)^2)*e_previous(2,i)';

end

if t<5% & single == 1
    G1_x = G1_x + .5*(rand(2,1)-[0.5;0.5]);
    G2_x = G2_x + .5*(rand(2,1)-[0.5;0.5]);
end



%% Transformed dynamics
trans1 = (((a1^2*exp(s(1))) - (2*a1*A1) + (A1^2 * exp (-s(1)))) / ( A1*a1^2 -a1*A1^2));
trans2 = (((a2^2*exp(s(2))) - (2*a2*A2) + (A2^2 * exp (-s(2)))) / ( A2*a2^2 -a2*A2^2));

F1_1 = trans1 * F1(1,:);
F1_2 = trans2 * F1(2,:);

G1_ = trans2 * (cos(2*x1)+2);
G2_ = trans2 * (sin(4*x1^2)+2);

F1_ = [F1_1; F1_2];
G1_ = [0; G1_];
G2_ = [0; G2_];


%% Terms in weight update laws that are evaluated along the system trajectory x

% Basis vector: [s1^2, s1s2,s2^2]
phi = [s(1)^2 s(1)*s(2) s(2)^2];
phi_p = [2*s(1) 0; s(2) s(1); 0 2*s(2)]; % Jacobian of basis vector


%% 控制策略

% rng(0,'twister');
% r = randi([1 6],1);

mu1a=-0.5.*(R1\G1_')*phi_p'*W_optimal(1,:)';
mu2a=-0.5.*(R2\G2_')*phi_p'*W_optimal_unsafe(2,:)';


a=[1 3 4 5 6];
b=[1 2 3 4 5];
chosen = r_sum(a);

sigma = softmax(-chosen');

Prob=sigma';

S=randsrc(1,1,[b;Prob]);

C = Prob(S);

mu1i = -0.5.*(R1\G1_')*phi_p'*W_stack(S*2+1,:)';
mu2i = -0.5.*(R2\G2_')*phi_p'*W_stack_unsafe(S*2+2,:)'; % 机器变为unsafe

delta1=0.1;
delta2=0.9;
delta3=0.7;

if C<delta1
    alpha = 0;
elseif C<delta2
    alpha = delta3/(delta2-delta1)*C;
else
    alpha = delta3;
end
alpha

% mu12i = (1-alpha)*mu1i + alpha*mu2i;
s_ = F1_+ G1_*(1-alpha)*mu1i + G2_*alpha*mu2c;

% s_ = F1_+ G1_*mu12i + G2_*mu12i;
x_ = [s_(1)/trans1; s_(2)/trans2];

%% Integrating function

zDot = [s_; % s
       WaH1_; % WaH
       WcH1_; % WcH
       WaH2_; % WaH
       WcH2_; % WcH
       x_; % x
       ];
e_prev = zeros(2,10);
e_prev(1,1) = delta1;
e_prev(2,1) = delta2;
e_prev(1,2:10) = e_previous(1,1:9);
e_prev(2,2:10) = e_previous(2,1:9);

sigma = zeros(3,10);
sigma(:,1) = omega;
sigma(:,2:10) = sigma_previous(:,1:9);
end
































