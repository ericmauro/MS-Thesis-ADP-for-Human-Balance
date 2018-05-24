function [x_data,t_data] = bbadp()
% Balance Control Dynamic Model
% model-based LQR control with feedforward
clear all; close all; clc;
global A Ab B Bb D E D0 K_star L Lstar X0 X2 X3 X4 K0 random1 random2 random3 Ldag K
Nrand = 5;
range_rand = 100;
random1 = range_rand*rand(Nrand,1)*1-0.5*range_rand;
random2 = range_rand*rand(Nrand,1)*1-0.5*range_rand;
random3 = range_rand*rand(Nrand,1)*1-0.5*range_rand;
%% load model
angle_type = 1; % 1 for angles to vertical, 2 for segment angles
[A0,B0] = balance_models(1,angle_type); % Initial model, no cart
[A,B,C,D,Ab,Bb,Db] = balance_models(2,angle_type);

R = eye(3);
%Q = 0.001*diag([1 15 5 20 20 20]);
%Q = 1000*diag([0.048 0.989 131.64 0 0 0]);
Q = diag([20^2,100^2,50^2,0,0,0]);
%Q = diag(deg2rad([75,135,120,0,0,0]).^2);
%f0_max = 1000;
E = -1.4;
D0 = 6;
%% K0 
K0 = lqr(A0,B0,Q,R);
L0 = [0 0 0]';
eig(A-B*K0);
%% model-based LQR solution
K_star = lqr(A,B,Q,R);
eig(A-B*K_star);
% K0 = K_star+0.1*rand(3,6);
% eig(A-B*K0)
%% output regulation solution
AA = [A B;
      C zeros(3)];
BB = [-D;
       zeros(3,1)];
XUstar = inv(AA)*BB;
Xstar = XUstar(1:6);
Ustar = XUstar(7:9);

Lstar = Ustar+K_star*Xstar;

%% construct vector X
Null_C = null(C);
X0=zeros(6,1);
X2=[0;0;0;1;0;0];
X3=[0;0;0;0;1;0];
X4=[0;0;0;0;0;1];
%% Get trajectory for healthy individual on solid ground
if angle_type == 1
    x0 = [0.1 -0.2 -0.1 0 0 0 0 0]';
else
    x0 = [-0.2 -0.3 0 0 0 0 0 0]';
end
[x_sim1,t_sim1] = balance_normal_sim(A0,B0,zeros(2,6),zeros(2,3),x0,K0);

figure
subplot(3,1,1)
plot(t_sim1,x_sim1(:,1:3),'Linewidth',2);
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Angle $\theta$ [rad]','Fontsize',16);set(ylab,'Interpreter','latex');
legend({'Ankle','Knee','Hip'},'Fontsize',12)
title('Balance on Solid Ground')
xlim([0 1])
subplot(3,1,2)
plot(t_sim1,x_sim1(:,4:6),'Linewidth',2);
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Anglular Velocity $\dot{\theta}$ [rad/s]','Fontsize',16);set(ylab,'Interpreter','latex');
legend({'Ankle','Knee','Hip'},'Fontsize',12)
xlim([0 1])
subplot(3,1,3)
plot(t_sim1,x_sim1(:,7),'Linewidth',2);
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Cart Position (m)','Fontsize',16);set(ylab,'Interpreter','latex');
xlim([0 1])
%% Get trajectory for healthy individual on cart with no learning
[x_sim2,t_sim2] = balance_normal_sim(A,B,Ab,Bb,x0,K0);

figure
subplot(3,1,1)
plot(t_sim2,x_sim2(:,1:3),'Linewidth',2);
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Angle $\theta$ [rad]','Fontsize',16);set(ylab,'Interpreter','latex');
legend({'Ankle','Knee','Hip'},'Fontsize',12)
title('Balance on Cart Before Learning, No Force')
xlim([0 1])
subplot(3,1,2)
plot(t_sim2,x_sim2(:,4:6),'Linewidth',2);
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Anglular Velocity $\dot{\theta}$ [rad/s]','Fontsize',16);set(ylab,'Interpreter','latex');
legend({'Ankle','Knee','Hip'},'Fontsize',12)
xlim([0 1])
subplot(3,1,3)
plot(t_sim2,x_sim2(:,7),'Linewidth',2);
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Cart Position (m)','Fontsize',16);set(ylab,'Interpreter','latex');
xlim([0 1])

%% learning phase
index=1;errorK=[];errorP=[]; errorL=[];
K_1 = zeros(3,6); P_1 = zeros(6,6);
X_traj=[];
t_traj=[];
u0 = zeros(3,1); % for initialization
uh0 = u0;
T = 1/30;
K = K0;
Ldag = L0;
index=0;temp=0;INDEX=[];
x_mse_save = [];
dx_mse_save = [];
num_trials = 10;
num_learn = 0;
num_delay = 0;
push_type = 1; % 1 for just pushes from behind, 2 for pattern of forward/backward pushes
for j = 1:num_trials
    j
%x_save = [];
x0 = [0 0 0 0 0 0]';
xb = [0 0]';

%f0 = 0.1*rand(1);
f0 = 0.1+0.1*rand(1);
w0 = f0;
if push_type == 2
    if (mod(j,6) == 3)|(mod(j,6) == 5)
        f0 = -f0;
    end
    if j>0
        w0 = f0;
    end
end

Dee0=[];Dee2=[];Dee3=[];Dee4=[];
EE0=[];EE2=[];EE3=[];EE4=[];
WE0=[];WE2=[];WE3=[];WE4=[];
VE0=[];VE2=[];VE3=[];VE4=[];
X_traj=[];
t_traj=[];
for i=1:num_delay
    i;
    [t,X]=ode45(@mysys_delay,[(i-1)*T,i*T],[x0;f0;w0;xb]);
    X_traj=[X_traj;X];
    t_traj=[t_traj;t];
    x0=X(end,1:6)';
    f0=X(end,7);
    w0=X(end,8);
    xb=X(end,9:10)';
end
for i=num_delay+1:150
    %x_save = [x_save x0'];
    i;
    [t,X]=ode45(@mysys,[(i-1)*T,i*T],[x0;f0;
            kron(x0',x0')';kron(x0,uh0);kron(x0',w0')';kron(x0',x0')';kron(x0,uh0);kron(x0',w0')';kron(x0',x0')';
            kron(x0,uh0);kron(x0',w0')';kron(x0',x0')';kron(x0,uh0);kron(x0',w0')';w0;xb]);
    
    %X_traj=[X_traj;X(:,:)];
    X_traj = [X_traj;X(:,1:7),X(:,248:250)];
    t_traj=[t_traj;t];
    
    E0=X(end,1:6)-X(end,248)*X0'; %x-Xi*f at time t1
    E0_pre=X(1,1:6)-X(1,248)*X0';%x-Xi*f at time t0

    E2=X(end,1:6)-X(end,248)*X2';
    E2_pre=X(1,1:6)-X(1,248)*X2';
    
    E3=X(end,1:6)-X(end,248)*X3';
    E3_pre=X(1,1:6)-X(1,248)*X3';
    
    E4=X(end,1:6)-X(end,248)*X4';
    E4_pre=X(1,1:6)-X(1,248)*X4';
    
    Dee0=[Dee0;kron(E0,E0)-kron(E0_pre,E0_pre)]; %delta_xbar_xbar
    Dee2=[Dee2;kron(E2,E2)-kron(E2_pre,E2_pre)];
    Dee3=[Dee3;kron(E3,E3)-kron(E3_pre,E3_pre)];    
    Dee4=[Dee4;kron(E4,E4)-kron(E4_pre,E4_pre)];
    
    EE0=[EE0;X(end,8:43)-X(1,8:43)]; %Gamma_xbar_xbar
    WE0=[WE0;X(end,44:61)-X(1,44:61)]; %Gamma_xbar_w
    VE0=[VE0;X(end,62:67)-X(1,62:67)]; %Gamma_xbar_v

    EE2=[EE2;X(end,68:103)-X(1,68:103)];
    WE2=[WE2;X(end,104:121)-X(1,104:121)];
    VE2=[VE2;X(end,122:127)-X(1,122:127)];
    
    EE3=[EE3;X(end,128:163)-X(1,128:163)];
    WE3=[WE3;X(end,164:181)-X(1,164:181)];
    VE3=[VE3;X(end,182:187)-X(1,182:187)];
    
    EE4=[EE4;X(end,188:223)-X(1,188:223)];
    WE4=[WE4;X(end,224:241)-X(1,224:241)];
    VE4=[VE4;X(end,242:247)-X(1,242:247)];
    
      
    x0=X(end,1:6)';
    f0=X(end,7);
    w0=X(end,248);
    xb=X(end,249:250)';
end
% figure(1);
% plot(t_traj,X_traj(:,1:6),'Linewidth',2);
% xlab = xlabel('Time [s]','Fontsize',20);set(xlab,'Interpreter','latex')


% if ((rank([EE0,WE0,VE0])<45)||(rank([EE2,WE2,VE2])<45)||(rank([EE3,WE3,VE3])<45)||(rank([EE4,WE4,VE4])<45))
% error('The rank conditon is not satisfied!')
% end
% rank([EE0,WE0])
% if ((rank([EE0,WE0])<39))
% error('The rank conditon is not satisfied!')
% end


% Dee0=Dee0(:,[1:6,8:12,15:18,22:24,29:30,36]); 
% Dee2=Dee2(:,[1:6,8:12,15:18,22:24,29:30,36]);  
% Dee3=Dee3(:,[1:6,8:12,15:18,22:24,29:30,36]); 
% Dee4=Dee4(:,[1:6,8:12,15:18,22:24,29:30,36]); 



%index=0;temp=0;INDEX=[];
%
%for j=1:30
    K_save{j} = K;
    L_save{j} = Ldag;
    QK = Q+K'*R*K;                    % Update the Qk matrix
    
    Theta = [Dee0,-2*EE0*kron(eye(6),K')-2*WE0, -2*VE0];  % Left-hand side of
    % the key equation
    Xi = -EE0*QK(:);                 % Right-hand side of the key equation
    pp = pinv(Theta)*Xi;             % Solve the equations in the LS sense
    P = reshape(pp(1:6*6), [6, 6]);  % Reconstruct the symmetric matrix
    P = (P + P');
    
%     BPv = pp(end-(6*3-1):end);
    BPv = pp(54-(6*3-1):54);
    K = inv(R)*reshape(BPv,3,6);% Get the improved gain matrix
errorK=[errorK;norm(K-K_star)];
errorP=[errorP;norm(P-P_1)];
INDEX=[INDEX;j];
%{
if (norm(K-K_1)<1e-6) 
    disp('meet the converngent requirement and the norm of the P_k-P_{k-1} is')
    K
    P 
%     barAX0

    break
end 
%}
K_1=K;
P_1=P;
% if j==50
%     error('The P and K are not convergent!')
% end
if j == 1
    X1 = X_traj;
    t1 = t_traj;
end

x_mse = [immse(X_traj(:,1),zeros(length(X_traj(:,1)),1)); ...
            immse(X_traj(:,2),zeros(length(X_traj(:,2)),1)); ...
            immse(X_traj(:,3),zeros(length(X_traj(:,3)),1))];
dx_mse = [immse(X_traj(:,4),zeros(length(X_traj(:,4)),1)); ...
            immse(X_traj(:,5),zeros(length(X_traj(:,5)),1)); ...
            immse(X_traj(:,6),zeros(length(X_traj(:,6)),1))];
%x_mse = immse(X_traj(:,1:3),zeros(size(X_traj(:,1:3))));
%dx_mse = immse(X_traj(:,4:6),zeros(size(X_traj(:,4:6))));
x_mse_save = [x_mse_save x_mse];
dx_mse_save = [dx_mse_save dx_mse];
%end
if j==1
B_estimated=inv(P)*K'; %B_estimated is constructed by online data

D_estimated = inv(P)*pp(55:60);

%% find barAX2
Xi = -EE2*QK(:);
Theta = [Dee2,-2*EE2*kron(eye(6),K')-2*WE2, -2*VE2];  % Left-hand side of the key equation
pp = pinv(Theta)*Xi;             % Solve the equations in the LS sense
barAX2 = D_estimated - inv(P)*pp(55:60);
AX2 = inv(P)*pp(55:60) - D_estimated ;
AX2_bar_real = -A*X2;

%% find barAX3
Xi = -EE3*QK(:);
Theta = [Dee3,-2*EE3*kron(eye(6),K')-2*WE3, -2*VE3];  % Left-hand side of the key equation
pp = pinv(Theta)*Xi;             % Solve the equations in the LS sense
barAX3 = D_estimated - inv(P)*pp(55:60);
AX3 = inv(P)*pp(55:60) - D_estimated ;
AX3_bar_real = -A*X3;

%% find barAX4
Xi = -EE4*QK(:);
Theta = [Dee4,-2*EE4*kron(eye(6),K')-2*WE4, -2*VE4];  % Left-hand side of the key equation
pp = pinv(Theta)*Xi;             % Solve the equations in the LS sense
barAX4 = D_estimated - inv(P)*pp(55:60);
AX4 = inv(P)*pp(55:60) - D_estimated ;
AX4_bar_real = -A*X4;

%% 
capA = [barAX2 barAX3 barAX4  -B_estimated];% eq (41) on Weinan's IFAC paper
capX = inv(capA)*(D_estimated);
Xdag = X0+capX(1)*X2+capX(2)*X3+capX(3)*X4;
Udag=capX(4:6);
alpha2 = capX(1);
alpha3 = capX(2);
alpha4 = capX(3);
Ldag = Udag+K*Xdag;
%elseif j<num_learn
  %  Ldag = [0 0 0]';
 %   Udag = [0 0 0]';
 %   Xdag = [0 0 0]';
end
errorL=[errorL;norm(Ldag-Lstar)];
end
%% show result
disp(['K_', num2str(INDEX(end)), '=']);
disp(K);
disp('K*=')
disp(K_star)

disp(['P_', num2str(INDEX(end)), '=']);
disp(P);

disp(['Xdag =']);
disp(Xdag);
disp('X*=')
disp(Xstar)

disp(['Udag =']);
disp(Udag);
disp('U*=')
disp(Ustar)

disp(['Ldag =']);
disp(Ldag)
disp('L*=')
disp(Lstar)
% learning ends
% show plots 1

figure
subplot(3,1,1)
hold on
plot(t1,X1(:,1:3),'Linewidth',2);
plot([1.5 1.5],ylim,'k-');
hold off
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Angle $\theta$ [rad]','Fontsize',16);set(ylab,'Interpreter','latex');
legend({'Ankle','Knee','Hip'},'Fontsize',12)
title('Trial 1 States')
xlim([0 5])
subplot(3,1,2)
hold on
plot(t1,X1(:,4:6),'Linewidth',2);
plot([1.5 1.5],ylim,'k-');
hold off
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Anglular Velocity $\dot{\theta}$ [rad/s]','Fontsize',16);set(ylab,'Interpreter','latex');
legend({'Ankle','Knee','Hip'},'Fontsize',12)
xlim([0 5])
subplot(3,1,3)
hold on
plot(t1,X1(:,7),'Linewidth',2);
plot([1.5 1.5],ylim,'k-');
hold off
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Applied Force (N)','Fontsize',16);set(ylab,'Interpreter','latex');
xlim([0 5])

figure
subplot(2,1,1)
%hold on
plot(t1,X1(:,7),'b','Linewidth',2);
%plot(t1,X1(:,248),'r','LineWidth',2);
%hold off
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Applied Force','Fontsize',16);set(ylab,'Interpreter','latex');
%legend('Actual force','Estimated force')
title('Trial 1, force and cart position')
xlim([0 5])
subplot(2,1,2)
%plot(t1,X1(:,249),'Linewidth',2);
plot(t1,X1(:,9),'LineWidth',2);
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Cart Position','Fontsize',16);set(ylab,'Interpreter','latex');
xlim([0 5])

figure
subplot(1,2,1);plot(INDEX,errorP,'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
xlab = xlabel('Number of Iterations','Fontsize',14);set(xlab,'Interpreter','latex');
leg = legend('$|P_{j}-P^{*}|$');set(leg,'FontSize',18);set(leg,'Interpreter','latex');
set(gca,'FontSize',16);
xlim([0 10])

subplot(1,2,2);plot(INDEX,errorK,'--rs','LineWidth',2, 'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
xlab = xlabel('Number of Iterations','Fontsize',14);set(xlab,'Interpreter','latex');
leg = legend('$|K_{j}-K^{*}|$');set(leg,'FontSize',18);set(leg,'Interpreter','latex');
set(gca,'FontSize',16);
xlim([0 10])

figure
plot(INDEX, errorL,'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
xlab = xlabel('Number of Iterations','Fontsize',14);set(xlab,'Interpreter','latex');
leg = legend('$|L_{j}-L^{*}|$');set(leg,'FontSize',18);set(leg,'Interpreter','latex');
set(gca,'FontSize',16);
xlim([0 num_trials])


figure
%subplot(3,1,1)
plot(INDEX, sum(x_mse_save,1)./3,'bo-','LineWidth',2)
xlabel('Trial #')
ylabel('Simulation MSE')
title('Angles')

figure
%subplot(3,1,1)
plot(INDEX, sum(dx_mse_save,1)./3,'ro-','LineWidth',2)
xlabel('Trial #')
ylabel('Simulation MSE')
title('Angular Velocities')


figure
subplot(2,1,1)
plot(t_traj,X_traj(:,1:3),'Linewidth',2);
xlabel('Time','Fontsize',16);
ylab = ylabel('Angle $\theta$ [rad]','Fontsize',16);set(ylab,'Interpreter','latex');
legend({'Ankle','Knee','Hip'},'Fontsize',12)
title(['Trial ',num2str(num_trials),' States'])
xlim([0 5])
subplot(2,1,2)
plot(t_traj,X_traj(:,4:6),'Linewidth',2);
xlabel('Time','Fontsize',16);
ylab = ylabel('Anglular Velocity $\dot{\theta}$ [rad/s]','Fontsize',16);set(ylab,'Interpreter','latex');
legend({'Ankle','Knee','Hip'},'Fontsize',12)
xlim([0 5])

figure
subplot(2,1,1)
%hold on
plot(t_traj,X_traj(:,7),'b','Linewidth',2);
%plot(t_traj,X_traj(:,248),'r','LineWidth',2);
%hold off
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Applied Force','Fontsize',16);set(ylab,'Interpreter','latex');
%legend('Actual force','Estimated force')
title(['Trial ',num2str(num_trials),' , force and cart position'])
xlim([0 5])
subplot(2,1,2)
%plot(t_traj,X_traj(:,249),'Linewidth',2);
plot(t_traj,X_traj(:,9),'Linewidth',2)
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Cart Position','Fontsize',16);set(ylab,'Interpreter','latex');
xlim([0 5])

%%

T = 1/30;
j = num_trials;
for j = 1:size(K_save,2)
    X_traj_app=[];
    t_traj_app=[];
    x0 = [0 0 0 0 0 0]';
    
    %f0 = 0.1+0.1*rand(1);
    f0 = 0.1;
    w0 = f0;
    xb = [0 0]';
    %load('balance_board_parameters','x0');
    K = K_save{j};
    L = L_save{j};
    for i=1:num_delay
        [t,X]=ode45(@mysys_delay,[(i-1)*T,i*T],[x0;f0;w0;xb]);
        X_traj_app=[X_traj_app;X];
        t_traj_app=[t_traj_app;t];
        x0=X(end,1:6)';
        f0=X(end,7);
        w0=X(end,8);
        xb=X(end,9:10)';
    end
    for i=num_delay+1:150
        [t,X]=ode45(@mysys_init,[(i-1)*T,i*T],[x0;f0;w0;xb]);
        X_traj_app=[X_traj_app;X];
        t_traj_app=[t_traj_app;t];
        x0=X(end,1:6)';
        f0=X(end,7);
        w0=X(end,8);
        xb=X(end,9:10)';
    end
   t_data{j} = t_traj_app;
   x_data{j} = X_traj_app(:,1:6);
end

figure
subplot(3,1,1)
hold on
plot(t_traj_app,X_traj_app(:,1:3),'Linewidth',2);
plot([1.5 1.5],ylim,'k-');
hold off
xlabel('Time','Fontsize',16);
ylab = ylabel('Angle $\theta$ [rad]','Fontsize',16);set(ylab,'Interpreter','latex');
legend({'Ankle','Knee','Hip'},'Fontsize',12)
title('Post-Learning Trial States')
xlim([0 5])
subplot(3,1,2)
hold on
plot(t_traj_app,X_traj_app(:,4:6),'Linewidth',2);
plot([1.5 1.5],ylim,'k-');
hold off
xlabel('Time','Fontsize',16);
ylab = ylabel('Anglular Velocity $\dot{\theta}$ [rad/s]','Fontsize',16);set(ylab,'Interpreter','latex');
legend({'Ankle','Knee','Hip'},'Fontsize',12)
xlim([0 5])
subplot(3,1,3)
hold on
plot(t_traj_app,X_traj_app(:,7),'Linewidth',2);
plot([1.5 1.5],ylim,'k-');
hold off
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Applied Force (N)','Fontsize',16);set(ylab,'Interpreter','latex');
xlim([0 5])

figure
subplot(2,1,1)
%hold on
plot(t_traj_app,X_traj_app(:,7),'b','Linewidth',2);
%plot(t_traj_app,X_traj_app(:,8),'r','LineWidth',2);
%hold off
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Applied Force','Fontsize',16);set(ylab,'Interpreter','latex');
%legend('Actual force','Estimated force')
title('Post-Learning Trial, force and cart position')
xlim([0 5])
subplot(2,1,2)
plot(t_traj_app,X_traj_app(:,9),'Linewidth',2);
xlab = xlabel('Time [s]','Fontsize',16);set(xlab,'Interpreter','latex');
ylab = ylabel('Cart Position','Fontsize',16);set(ylab,'Interpreter','latex');
xlim([0 5])
%{
figure(4);
subplot(2,1,1)
plot([t_traj;t_traj(end)+t_traj_app],[X_traj(:,1:3);X_traj_app(:,1:3)],'Linewidth',2);
xlabel('Time','Fontsize',16);
ylab = ylabel('Angle $\theta$ [rad]','Fontsize',16);set(ylab,'Interpreter','latex');
legend('Ankle','Knee','Hip')
subplot(2,1,2)
plot([t_traj;t_traj(end)+t_traj_app],[X_traj(:,4:6);X_traj_app(:,4:6)],'Linewidth',2);
xlabel('Time','Fontsize',16);
ylab = ylabel('Anglular Velocity $\dot{\theta}$ [rad/s]','Fontsize',16);set(ylab,'Interpreter','latex');
legend('Ankle','Knee','Hip')
%}
end

% System subfunction
function dX=mysys(t,X) 
global A Ab B Bb D E D0 X0 X2 X3 X4 K_star Lstar K0 K random1 random2 random3 Ldag
x=X(1:6);
f=X(7);
w=X(248);
xb=X(249:250);

e0=x-X0*w;
% e1=x-X1*f;
e2=x-X2*w;
e3=x-X3*w;
e4=x-X4*w;


% u=sum(0.2*sin(random*t));
% u=sum(0.05*sin(random*t));
% u = - K_star * x + Lstar;

u = -K*x+Ldag*w+[sum(0.01*sin(random1*t)); sum(0.01*sin(random2*t)); sum(0.01*sin(random3*t))];
uh = u-Ldag*w;

dx=A*x+B*u+D*f;
df=-40*(t-2)*f*(t<2)*(t>1.5) + E*(f-D0)*(t>=2);
%df=E*(f-D0);
%dw=E*(w-D0);
%df = 0;
dw = df;
dxb=Ab*x+Bb*u+[0 1; 0 0]*xb;

de0e0=kron(e0',e0')'; %36
dwe0=kron(e0',uh')'; %18
dve0=kron(e0',w')'; %6

de2e2=kron(e2',e2')';
dwe2=kron(e2',uh')';
dve2=kron(e2',w')';

de3e3=kron(e3',e3')';
dwe3=kron(e3',uh')';
dve3=kron(e3',w')';

de4e4=kron(e4',e4')';
dwe4=kron(e4',uh')';
dve4=kron(e4',w')';


dX=[dx;
    df;
    de0e0; 
    dwe0; 
    dve0;    
    de2e2 
    dwe2; 
    dve2;
    de3e3 
    dwe3; 
    dve3;    
    de4e4
    dwe4; 
    dve4; 
    dw;
    dxb;
    ];
end


function dX=mysys_init(t,X)
% global A B C F K E Xdag Udag X0 X1 X2
global A Ab B Bb D E D0 K L Ldag Lstar K_star
x=X(1:6);
f=X(7);
w=X(8);
xb=X(9:10);
% u = zeros(3,1);
%u = - K_star * x + Lstar*w;
u = - K*x + L*w;
dx=A*x+B*u+D*f;
dxb=Ab*x+Bb*u+[0 1; 0 0]*xb;
df=-40*(t-2)*f*(t<2)*(t>1.5) + E*(f-D0)*(t>=2);
%df=E*(f-D0);
%dw=E*(w-D0);
%df = 0;
dw = df;
dX=[dx;df;dw;dxb];
end

function dX=mysys_delay(t,X)
% global A B C F K E Xdag Udag X0 X1 X2
global A Ab B Bb D E D0 Ki L Ldag Lstar K_star
x=X(1:6);
f=X(7);
w=X(8);
xb=X(9:10);
%u = zeros(3,1);
u = Ldag*w;
%u = - K_star * x + Lstar*w;
%u = - K * x + L*w;
dx=A*x+B*u+D*f;
dxb=Ab*x+Bb*u+[0 1; 0 0]*xb;
df=-40*(t-2)*f*(t<2)*(t>1.5) + E*(f-D0)*(t>=2);
%df=E*(f-D0);
%dw=E*(w-D0);
dw = df;
dX=[dx;df;dw;dxb];
end