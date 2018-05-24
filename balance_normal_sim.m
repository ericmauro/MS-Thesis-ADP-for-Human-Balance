function [X_traj,t_traj] = balance_normal_sim(A,B,Ab,Bb,x0,K)
% Eric Mauro
% Balance Control Dynamic Model
% Apr-10-2018
% 
% Simulate human balance as inverted pendulum on a cart.
% Inputs
% ---------
%   A:  pendulum system matrix
%   B:  pendulum input matrix
%   Ab: board (cart) system matrix
%   Bb: board (cart) input matrix
%   x0: initial state conditions
%   K:  feedback matrix
%
% Outputs:
% ---------
%   X_traj: state trajectories
%   t_traj: corresponding trajectory time instances

X_traj=[];
t_traj=[];
u0 = zeros(3,1);

T = 1/30;   % Sampling period
for i=1:150
        
        [t,X]=ode45(@mysys,[(i-1)*T,i*T],x0);
        X_traj=[X_traj;X];
        t_traj=[t_traj;t];
        x0=X(end,:)';
end
function dX=mysys(t,X)
x=X(1:6);
xb=X(7:8);
u = - K * x;
dx=A*x+B*u;
dxb=Ab*x+Bb*u+[0 1; 0 0]*xb;
dX=[dx;dxb];
end
end

