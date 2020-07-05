%% Periodogram of final SDE for different frequencies
%  dx = (x - x^3 - V1*sin(omega*t))*dt + kappa*dW

% Time data
Fs = 1000; % Sampling frequency, 1kHz
dt = 1/Fs; % Sampling period
T = 2000; % Final time (in seconds)
t = 0:dt:T-dt; % Time vector
L = length(t); % Length of the signal; # of samples

% Particle position data
x = zeros(1,L); % Particle's position
x(1) = 0.5;   % Initial condition

% Potential data
V0 = @(y) -y.^2./2 + y.^4./4; % Potential
dV0 = @(y) -y + y^3;
DV0 = V0(0) - V0(1); % Barrier
V1 =  DV0*0.2; % Small compared to barrier
F = [1/200,1/100,1/50]; % Frequency of the modulated potential
df1 = @(s) V1*sin(2*pi*F*s); % Modulated potential

kappa = 0.5; % Noise strength

M = 50; % Number of realisations
xM = zeros(M,L);

for k = 1:3
    for m = 1:M
        % EULER MARUYAMA method
        for n = 1:L-1
        % POTENTIAL PERIODIC NOISE
            x(n+1) = x(n) - dt*( dV0(x(n)) + V1*sin(2*pi*F(k)*t(n+1)) ) + sqrt(kappa*dt)*randn; 
        end 
        xM(m,:) = x;
    end
    xmean = mean(xM,1);
    subplot(3,1,k)
        periodogram(xmean,rectwin(L),L,Fs)
        xlim([0,.1])
        title(['$\omega_{s} = $ ', num2str(F(k))],'Interpreter', 'latex','FontSize',15)
end
sgtitle('Periodograms of the final SDE for different $\omega_{s}$','Interpreter', 'latex','FontSize',20)