%% Periodogram of final SDE for different kappas
%  dx = (x - x^3 - V1*sin(omega*t))*dt + kappa*dW

% Time data
Fs = 1000; % Sampling frequency, 1kHz
dt = 1/Fs; % Sampling period
T = 2000; % Final time (in seconds)
t = 0:dt:T-dt; % Time vector
L = length(t); % Length of the signal; # of samples

% Particle position data
x = zeros(1,L); % Particle's position
x(1) = .5; % Initial condition

% Potential data
V0 = @(y) -y.^2./2 + y.^4./4; % Potential
dV0 = @(y) -y + y^3;
DV0 = V0(0) - V0(1); % Barrier
V1 =  DV0*0.2; % Small compared to barrier
F = 1/100; % Frequency of the modulated potential
omega = 2*pi*F; % Angular frequency
df1 = @(s) V1*sin(omega*s); % Modulated potential

% Noise strength
k1 = 2*(DV0-V1)/log(1/(omega*sqrt(2))); %k1
k2 = 2*(DV0+V1)/log(1/(omega*sqrt(2))); %k2
kappa = [k1,(k1+k2)/2,k2]; 
Lk = length(kappa);

M = 50; % Number of realisations
xM = zeros(M,L); % Mean of x

for k = 1:Lk
    for m = 1:M
        % EULER MARUYAMA method to get x(t)
        for n = 1:L-1
        % POTENTIAL PERIODIC NOISE
            x(n+1) = x(n) - dt*( dV0(x(n)) + df1(t(n+1)) ) + sqrt(kappa(k)*dt)*randn; 
        end 
        xM(m,:) = x;
    end
    xmean = mean(xM,1);
    subplot(Lk,2,2*k-1)
        plot(t,x,t,xmean)
        title(['$\kappa = $ ', num2str(kappa(k))],'Interpreter', 'latex')
        ylabel('Position')
        xlabel('Time (s)')
        legend('x(t)',['Mean of ',num2str(M),' realisations'])
    subplot(Lk,2,2*k)
        periodogram(xmean,rectwin(L),L,Fs)
        title('Periodogram')
        xlim([0,.05])
end
sgtitle('Position of the final SDE and periodogram for different $\kappa$','Interpreter', 'latex')