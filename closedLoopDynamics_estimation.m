function [zDot, e_prev,sigma] = closedLoopDynamics_estimation(t,z)
global e_previous sigma_previous j level W_stack W_optimal
% Control Penalty Matrix  
R1=0.1*10;
R2=R1/2;

a1 = -1.3; A1 = 0.5; %Barrier Boundaries for sii(1), sii(1)-coordinate values will be remained within (a1,A1) 
a2 = -3.1; A2 = 0.5; %Barrier Boundaries for sii(2), sii(2)-coordinate values will be remained within (a2,A2) 

s = z(1:2,1); % row 1-2 in z vector


%% Dynamics

%%transforming initial s-coordinates to initial x-coordinates. 
x1 = a1*A1*((exp(s(1)))-1)/((a1*exp(s(1)))-A1);
x2 = a2*A2*((exp(s(2)))-1)/((a2*exp(s(2)))-A2);


%% Original dynamics without the parameters
F1 = [x2; -x2 - 0.5*x1 + 0.25*x2*(cos(2*x1)+2)^2 + 0.25*x2*(sin(2*x1)+2)^2];
G1 = [0; (cos(2*x1)+2)];
G2 = [0; (sin(4*x1^2)+2)];

%%Transformed dynamics
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

rng(0,'twister');
r = randi([1 6],1);

mu1a=-0.5.*(R1\G1_')*phi_p'*W_stack(r*2+1,:)';
mu2a=-0.5.*(R2\G2_')*phi_p'*W_stack(level*2+2,:)';

s_ = F1_+ G1_*mu1a + G2_*mu2a;
x_ = [s_(1)/trans1; s_(2)/trans2];


%% Integrating function

zDot = [s_; % s
       x_ % x
       ];

end
































