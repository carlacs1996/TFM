%% 2.3 WHITE NOISE
%  dx = (x - x^3)*dt + sqrt(kappa)*dW
clc
clear

% Setting initial data
dt = 0.01;
T = 2000;
t = 0:dt:T;
N = T/dt;
x = zeros(3,N+1);
x(:,1) = 1;

% Potential data
V0 = @(y) -y.^2./2 + y.^4./4; % Potential
dV0 = @(y) -y + y^3;
ddV0 = @(y) -1 + 3*y^2;
DV0 = V0(0) - V0(1); % Barrier

% Noise strength
kappa = [0.2; 0.5; 0.8];
 
% EULER MARUYAMA 
tau = zeros(3,1); % Expected residence time
for k = 1:3
    for n = 1:N
        x(k,n+1) = x(k,n) - dt*dV0(x(k,n)) + sqrt(kappa(k)*dt)*randn;
    end
    tau(k) = 2*pi*exp(2*DV0/kappa(k))/sqrt( abs(ddV0(1)*ddV0(0)) );

    subplot(3,1,k)
    plot(t,x(k,:),'Color','#EDB120')
    title(['$\kappa = $ ', num2str(kappa(k)), '; average residence = ', num2str(round(tau(k))),' s'],'Interpreter', 'latex')
    xlabel('Time (s)','Interpreter', 'latex')
    ylabel('Position $x(t)$','Interpreter', 'latex')
    xlim([0,1000])
    grid on
    set(gca,'FontSize',20)
end
sgtitle('Potential and white noise','Interpreter', 'latex','FontSize',20)
%%
positive = x>=0; % Positive values of x
sp = sum(positive,2); % Number of positive values of x
Esp = sp./(N+1); % Average of values of x being positive
