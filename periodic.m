%% POTENTIAL AND PERIODIC FORCE
%  dx = (x - x^3 - V1*sin(omega*t))*dt

% Setting initial data
dt = 0.01; % Sampling interval
T = 2000; % Final time (seconds)
t = 0:dt:T;
N = T/dt;
x = zeros(1,N+1); % Particle's position
x(1) = 0.1; 

% Potential data
V0 = @(y) -y.^2./2 + y.^4./4; % Potential
dV0 = @(y) -y + y^3;
DV0 = V0(0) - V0(1); % Barrier height

% Periodic data
V1 =  DV0*0.2; % Small compared to barrier
F = 1/100; % Frequency of the modulated periodic force
omega = 2*pi*F; % Angular frequency
df1 = @(s) V1*sin(omega*s);

% EULER method
for n = 1:N
    x(n+1) = x(n) - ( dV0(x(n)) + df1((n-1)*dt) )*dt;
end

plot(t,x,'Color','#EDB120','LineWidth',3)
    title('Modulated periodic potential','Interpreter', 'latex')
    xlabel('Time (s)','Interpreter', 'latex')
    ylabel('Position $x(t)$','Interpreter', 'latex')
    xlim([0,500])
    grid on
    set(gca,'FontSize',20)