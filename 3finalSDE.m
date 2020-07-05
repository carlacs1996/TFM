%% 2.4 FINAL SDE
%  dx = (x - x^3 - V1*sin(omega*t))*dt + sqrt(kappa)*dW
clc
clear

% Initial data
dt = 0.01; % Sampling interval
T = 2000; % Final time (seconds)
t = 0:dt:T;
N = T/dt;
x = zeros(3,N+1);
x(:,1) = 1; 

%Potential data
V0 = @(y) -y.^2./2 + y.^4./4; %Potential
dV0 = @(y) -y + y^3;
DV0 = V0(0) - V0(1); % Barrier

% Periodic data
V1 =  DV0*0.3; % Small compared to barrier
F = 1/100; % Frequency of the modulated periodic force
omega = 2*pi*F;
f1 = @(y,s) V1*y*sin(omega*s);
df1 = @(s) V1*sin(omega*s);

% Noise strength
kappa = [0.2; 0.5; 0.8];

% EULER MARUYAMA method
for k = 1:3
    for n = 1:N
        x(k,n+1) = x(k,n) - dt*( dV0(x(k,n)) + df1(t(n)) ) + sqrt(kappa(k)*dt)*randn;
    end
    
    subplot(3,1,k)
        plot(t,x(k,:),'Color','#EDB120')
        title(['$\kappa = $', num2str(kappa(k))],'Interpreter', 'latex')
        xlabel('Time','Interpreter', 'latex')
        ylabel('Position $x(t)$','Interpreter', 'latex')
        grid on
        xlim([0,1000])
        set(gca,'FontSize',20)
end
sgtitle('Periodic modulated potential with white noise','Interpreter', 'latex','FontSize',20)