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

s_ = F1_x+ G1_x*mu1c + G2_x*mu2c;

%% Integrating function

zDot = [s_; % s
       WaH1_; % WaH
       WcH1_; % WcH
       WaH2_; % WaH
       WcH2_; % WcH
       ];

%% 下面的是存储历史误差值和σ
e_prev = zeros(2,10);
e_prev(1,1) = delta1x;
e_prev(2,1) = delta2x;
e_prev(1,2:10) = e_previous(1,1:9);
e_prev(2,2:10) = e_previous(2,1:9);

sigma = zeros(3,10);
sigma(:,1) = omega;
sigma(:,2:10) = sigma_previous(:,1:9);
end






