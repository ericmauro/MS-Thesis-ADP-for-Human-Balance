function [A,B,C,D,Ab,Bb,Db] = balance_models(num_case,angle_type)
% Eric Mauro
% Balance Control Dynamic Model
% Apr-09-2018
% -----------------------------
% A, B, C, D are system matrices in 
%   \dot{x} = Ax + Bu + Df
%         y = Cx
% There are 2 cases that affect the dynamics:
%   1. Stationary, no cart
%   2. w/ cart
% There are also 2 types of models:
%   1. Angles referenced to vertical axis
%   2. Segment angles referenced to body

%% Physical Parameters
g = 9.81; % Gravitational constant (m/s/s)

m = [4;          % Mass of calf/lower leg (kg)
    7;           % Mass of thigh (kg)
    48];         % Mass of torso + head (kg)

L = [0.6;        % Length of calf/lower leg (m)
    0.5;         % Length of thigh (m)
    0.8];        % Length of torso + head (m)

c = [0.567;      % Ratio of l1, distance from ankle to CoM for calf
    0.567;       % Ratio of l2, distance from knee to CoM for thigh
    0.4829];     % Ratio of l3, distance to CoM of torso + head

I = [0.12;       % Moment of Inertia for calf (kg*m^2)
    0.1458;      % Moment of Inertia for thigh (kg*m^2)
    2.26];       % Moment of Inertia for torso + head (kg*m^2)

m0 = 2;          % Mass of board 
m_sum = m0+sum(m);

%% Linearized Model
% Linearized around upright position (theta = 0)
% cos -> 1, sin -> theta, \dot{theta}^2 -> 0
h = [m(1)*c(1)^2+m(2)+m(3); m(2)*c(2)^2+m(3); m(3)*c(3)^2]; % Constant
k = [m(1)*c(1)+m(2)+m(3); m(2)*c(2)+m(3); m(3)*c(3)];       % Constant

M = [I(1)+h(1)*L(1)^2 L(1)*L(2)*k(2) L(1)*L(3)*k(3);
    L(1)*L(2)*k(2) I(2)+h(2)*L(2)^2 L(2)*L(3)*k(3);
    L(1)*L(3)*k(3) L(2)*L(3)*k(3) I(3)+h(3)*L(3)^2]; % \ddot{theta} coeff. matrix

G = g.*[L(1)*k(1) 0 0; 
        0 L(2)*k(2) 0; 
        0 0 L(3)*k(3)]; % theta coeff. matrix
    
W = k.*L;

if angle_type == 1
    F = eye(3);
elseif angle_type == 2
    F = [1 0 0; 1 1 0; 1 1 1];
end

%% Set up State-Space Model and Simulation
if num_case == 1
    A = [zeros(3) eye(3); 
    inv(F)*inv(M)*G*F zeros(3)];
    B = [zeros(3); inv(F)*inv(M)];
    C = [eye(3) zeros(3,3)]; % Output matrix
    D = zeros(6,1); % Disturbance
elseif num_case == 2
    A = [zeros(3) eye(3); 
    inv(F)*inv(M-(1/m_sum)*W*W')*G*F zeros(3)];
    B = [zeros(3); inv(F)*inv(M-(1/m_sum)*W*W')];
    C = [eye(3) zeros(3,3)]; % Output matrix
    D = [zeros(3,1); inv(F)*inv(M-(1/m_sum)*W*W')*((-1/m_sum)*W)]; % Disturbance
    Ab = [zeros(1,6); W'*inv(M)*G*F/(m_sum+W'*inv(M)*W) zeros(1,3)];
    Bb = [zeros(1,3); W'*inv(M)/(m_sum+W'*inv(M)*W)];
    Db = 1/(m_sum+W'*inv(M)*W);
else
    disp('Unrecognized input')
end

end

