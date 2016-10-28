clear all;
% ---- Define System
v=0.8;
l=0.8;
m=1.5;
h=0.01;
R=10;
A=[0 1;(v*h/(m*(l^2)))-1 2-(v*h/(m*(l^2)))];
B=[0;(h^2)/(m*l^2)];
C=[0 1];

% ---- Trajectory
t_st=(0:0.01:10);
l=length(t_st);
for i=0:l-1
r(i+1)=t_st(i+1)^3*(4-0.3*t_st(i+1))*0.01;
end
figure(1)
plot(t_st,r);title('refrence')

% ---- Solve the riccati equation backwards
K(:,:,1001)=zeros(2,2); % K(N)
G(1001,:)=zeros(1,2); % Gain
for t=1:1000
    K(:,:,1001-t)=A'*K(:,:,1001-t+1)*A + C'*C - A'*K(:,:,1001-t+1)*B* ... 
        inv(B'*K(:,:,1001-t+1)*B + R)*B'*K(:,:,1001-t+1)*A;
    G(1001-t,:) = inv(B'*K(:,:,1001-t)*B+R)*B'*K(:,:,1001-t)*A; % Gain
end

% ---- Intialize
e=zeros(1,1001); % initial error
u=zeros(1,1000); % intialize input for previous trial
u(1)=0.6;
x(:,1)=zeros(2,1); % initial current position
u_k(1,1)=0.6; % intial current input
x1(:,1001)=zeros(2,1); % initialize states for previous trial

% ---- Run for 10 iterations
for k=1:10
    % ---- Solving eta backwards
    eta(:,1001)=zeros(2,1); % eta(N)
    for m=1:1000
        eta(:,1001-m) = inv(eye(2,2) + (K(:,:,1001-m) * B * (1/R) * B'))* ...
            (A'* eta(:,1001-m+1) + C' * e(1001-m+1));
    end
    for t=1:1000
        % ---- Update
        u_k(k,t) = u(t) - (G(t,:) * (x(:,t)-x1(:,t))) + ((1/R) * B' * eta(:,t)); % Input
        x(:,t+1) = A*x(:,t) + B*u_k(k,t); % States
        y(k,t)=C*x(:,t); % Output
        e(t)=r(t)-y(k,t); % Store Error
        x1(:,t)=x(:,t); % Store states for next trial
        u(t)=u_k(k,t); % Store input for next trial
    end
    x1(:,t+1)=x(:,t+1);
end
%% plot
figure(2);
subplot(2,1,1)% ---- output
plot(t_st(1:end-1), y(2,:), t_st(1:end-1), y(3,:), t_st(1:end-1), y(4,:), t_st(1:end-1), y(6,:), t_st(1:end-1), y(8,:), t_st(1:end-1), y(9,:));title('output for different iteration from 2-9 (rho=1)');legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 9');
subplot(2,1,2) % ---- input
plot(t_st(1:end-1), u_k(2,:), t_st(1:end-1), u_k(3,:), t_st(1:end-1), u_k(4,:), t_st(1:end-1), u_k(6,:), t_st(1:end-1), u_k(8,:), t_st(1:end-1), u_k(10,:));title('input for different iterations from 2-10 (rho=1)');legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 10')
figure(3)
subplot(2,1,1) % ---- output iteration 10 vs reference
plot(t_st(1:end-1), y(10,:), t_st(1:end-1), r(1:end-1));title('refernce v/s output iteration 10 (rho=10)');legend('output iteration 10','reference')
subplot(2,1,2) % ---- error
plot(t_st, e);title('error')
