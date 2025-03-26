function [zDot, e_prev,sigma] = unsafe(t,z)
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

% Optimizing Matrix/Integraing matrix >> This z matrix is the integration
% of the zDot matrix given at the end of this coding.  
s = z(1:2,1); % row 1-2 in z vector
s1 = s(1);
s2 = s(2);
WaH1 = z(3:5); % row 35-37 in z vector
WcH1 = z(6:8); % row 38-40 in z vector
WaH2 = z(9:11); % row 35-37 in z vector
WcH2 = z(12:14); % row 38-40 in z vector

%% Original dynamics without the parameters
F1 = [s2; -s2 - 0.5*s1 + 0.25*s2*(cos(2*s1)+2)^2 + 0.25*s2*(sin(2*s1)+2)^2];
G1 = [0; (cos(2*s1)-2)];
G2 = [0; (sin(4*s1^2)-2)];
%Known Parameters of the dynamics

F1_1 = F1(1,:);
F1_2 = F1(2,:);

G1_ = (cos(2*s1)+2);
G2_ = (sin(4*s1^2)+2);

F1_ = [F1_1; F1_2];
G1_ = [0; G1_];
G2_ = [0; G2_];

%% Terms in weight update laws that are evaluated along the system trajectory x

% Basis vector: [s1^2, s1s2,s2^2]
phi = [s(1)^2 s(1)*s(2) s(2)^2];
phi_p = [2*s(1) 0; s(2) s(1); 0 2*s(2)]; % Jacobian of basis vector
%phi_p_ = [2*s1 0; s2 s1; 0 2*s2]; % Jacobian of basis vector

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

%% Integrating function

zDot = [s_; % s
       WaH1_; % WaH
       WcH1_; % WcH
       WaH2_; % WaH
       WcH2_; % WcH
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







