function [zDot, e_prev,sigma] = closedLoopDynamics(t,z)
global e_previous sigma_previous j level W_stack W_optimal
% Control Penalty Matrix  
R1=0.1*10;
R2=R1/2;
% State Penalty Matrix  
Q=[10 0; 0 10]/10; 
% Tuning Gains
etac1 = 1;
etac2 = 1;
etaa1 = 30;
etaa2 = 30;
beta = 3;
v = 1;

single = 1;

M= 10;

a1 = -1.3; A1 = 0.5; %Barrier Boundaries for sii(1), sii(1)-coordinate values will be remained within (a1,A1) 
a2 = -3.1; A2 = 0.5; %Barrier Boundaries for sii(2), sii(2)-coordinate values will be remained within (a2,A2) 

% Optimizing Matrix/Integraing matrix >> This z matrix is the integration
% of the zDot matrix given at the end of this coding.  
s = z(1:2,1); % row 1-2 in z vector
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
%Known Parameters of the dynamics


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
%phi_p_ = [2*x1 0; x2 x1; 0 2*x2]; % Jacobian of basis vector

%% 更新律
%controller of the transformed system

mu1c = -0.5.*(R1\G1_')*phi_p'*WcH1;
mu2c = -0.5.*(R2\G2_')*phi_p'*WcH2;

if single == 0
    mu1a = -0.5.*(R1\G1_')*phi_p'*WaH1;
    mu2a = -0.5.*(R2\G2_')*phi_p'*WaH2;
elseif single == 1
    mu1a = -0.5.*(R1\G1_')*phi_p'*WcH1;
    mu2a = -0.5.*(R2\G2_')*phi_p'*WcH2;
else
    mu1a = WaH1'*phi';
    mu2a = WaH2'*phi';
end

%%
t0=30;
tao = 1;
if j == 0
    mu2a=-0.5.*(R2\G2_')*phi_p'*W_stack(level*2+2,:)';
%     mu1a = 0;
    if t < t0
        mu2a = mu2a*(1-exp(-tao*t));
    end
end
if j == 1
    mu1a=-0.5.*(R1\G1_')*phi_p'*W_stack(level*2+3,:)';
%     mu2a = 0;
    if t < t0
        mu1a = mu1a*(1-exp(-tao*t));
    end
end



%% Cost Function in the s-coordinate
r1 = s'*Q*s + mu1a'*R1*mu1a + mu2a'*R2*mu2a;
r2 = r1/2;

%% Definition
omega = phi_p*(F1_+G1_*mu1a+G2_*mu2a);
rho = (1+v*(omega'*omega))^2;

delta1 = WcH1'*omega + r1;
delta2 = WcH2'*omega + r2;

WcH1_ = -etac1*omega*delta1/rho;% + clWc1; % WcH
WcH2_ = -etac2*omega*delta2/rho;% + clWc2; % WcH

if single == 0 
    WaH1_ = - beta*(WaH1 - WcH1);
    WaH2_ = - beta*(WaH2 - WcH2);
elseif single == 1
    WaH1_ = - beta*(WaH1 - WcH1);
    WaH2_ = - beta*(WaH2 - WcH2);
else
    WaH1_ = -etaa1*phi'*(mu1a - mu1c);% + clWa1; % WaH
    WaH2_ = -etaa2*phi'*(mu2a - mu2c);% + clWa1; % WaH 
end


%% BE Extrapolation - Terms in weight update laws that are evaluated along arbitrarily selected trajectories sii

for i=1:M
    WcH1_ = WcH1_ - etac1*(sigma_previous(:,i)./(sigma_previous(:,i)'*sigma_previous(:,i)+1)^2)*e_previous(1,i)';
    WcH2_ = WcH2_ - etac2*(sigma_previous(:,i)./(sigma_previous(:,i)'*sigma_previous(:,i)+1)^2)*e_previous(2,i)';

end

if t<5% & single == 1
    G1_ = G1_ + .5*(rand(2,1)-[0.5;0.5]);
    G2_ = G2_ + .5*(rand(2,1)-[0.5;0.5]);
end

s_ = F1_+ G1_*mu1a + G2_*mu2a;
x_ = [s_(1)/trans1; s_(2)/trans2];

%% Integrating function

zDot = [s_; % s
       WaH1_; % WaH
       WcH1_; % WcH
       WaH2_; % WaH
       WcH2_; % WcH
       x_ % x
       ];

%% 下面的是存储历史误差值和σ
e_prev = zeros(2,10);
e_prev(1,1) = delta1;
e_prev(2,1) = delta2;
e_prev(1,2:10) = e_previous(1,1:9);
e_prev(2,2:10) = e_previous(2,1:9);

sigma = zeros(3,10);
sigma(:,1) = omega;
sigma(:,2:10) = sigma_previous(:,1:9);
end







