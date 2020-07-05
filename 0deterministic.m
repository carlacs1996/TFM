%% 2.1 POTENTIAL FORCE
%  dx = (x - x^3)*dt
clc
clear

% Setting initial data
dt = 0.01;
T = 20;
t = 0:dt:T;
N = T/dt;
x = zeros(1,N+1);
x(1) = 0.1; 

% Potential data
V0 = @(y) -y.^2./2 + y.^4./4; % Potential
dV0 = @(y) -y + y^3;
DV0 = V0(0) - V0(1); % Barrier length

% EULER method
for n = 1:N
    x(n+1) = x(n) - dt*dV0(x(n));
end

plot(t,x,'Color','#EDB120','LineWidth',3)
    title('Purely determinsitic potential','Interpreter', 'latex')
    xlabel('Time (s)','Interpreter', 'latex')
    ylabel('Position $x(t)$','Interpreter', 'latex')
    xlim([0,10])
    grid on
    set(gca,'FontSize',20)