%% ONLY DETERMINISTIC FORCE: SYMMETRIC DOUBLE WELL POTENTIAL
%  dx = (x - x^3)*dt

% Setting initial data
dt = 0.01; % Sampling interval
T = 20; % Final time (seconds)
t = 0:dt:T;
N = T/dt;
x = zeros(1,N+1); % Particle's position
x(1) = 0.1; 

% Potential data
dV0 = @(y) -y + y^3;

% EULER method
for n = 1:N
    x(n+1) = x(n) - dV0(x(n))*dt;
end

plot(t,x,'Color','#EDB120','LineWidth',3)
    title('Symmetric double well potential','Interpreter', 'latex')
    xlabel('Time (s)','Interpreter', 'latex')
    ylabel('Position $x(t)$','Interpreter', 'latex')
    xlim([0,10])
    grid on
    set(gca,'FontSize',20)