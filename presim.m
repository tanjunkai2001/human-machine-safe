%% Paper title "Safe Model-Based Reinforcement Learning for Systems with Parametric Uncertainties" 
%% arXiv:2007.12666
%% Coded by S M Nahid Mahmud, MS Grad Student, Oklahoma State University.
%% nahid.mahmud@okstate.edu
%% Two-state dynamical system (Section 7.1 of the journal paper)

clear
close all
clc
tic
%% PreSim
tic
a1 = -1.3; A1 = 0.5; %Barrier Boundaries for x1, x1-coordinate values will be remained within (a1,A1) 
a2 = -3.1; A2 = 0.5; %Barrier Boundaries for x2: x2-coordinate values will be remained within (a2,A2) 

P.n = 2; %number of states
P.m = 1; %number of dimension
P.l = 4; %number of unknown parameters
% P.x0 = [-1;-3]; % Initial x-coordinate values, it has to be within (a,A).
P.x0 = [-0.5;-2]; % Initial x-coordinate values, it has to be within (a,A).
%% 

%transforming initial x-coordinates to initial s-coordinates. 
s10 = (log ((A1/a1)*((a1-P.x0(1))/(A1-P.x0(1)))));
s20 = (log ((A2/a2)*((a2-P.x0(2))/(A2-P.x0(2)))));

P.s0 = [s10;s20]; % Initial s-coordinate values
P.W = 3; % number of the Actor Weight vector elements as wellas Critic Weight vector elements.
P.G =3;  % number of the Gamma vector elements. 

%Initial values for s = 1-2, WaH = 35:37, WcH = 38:40, x = 50:51, Int_G = 52:53
% WaH denotes Estimated Actor weight, WcH denotes Estimated Critic weight,
% Int_G denotes the integration of G function. 
% P.z0 = [P.s0;0.5*ones(4*P.W,1);P.x0];
P.z0 = [P.s0;ones(4*P.W,1)-0.5;P.x0];


e_prev = zeros(2,10);
sigma = zeros(3,10);