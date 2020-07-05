%% Single and average and periodogram for both noise and final SDE
%  dx = (x - x^3 - V1*sin(omega*t))*dt + sqrt(kappa)*dW

% Time data
Fs = 1000; % Sampling frequency, 1kHz
dt = 1/Fs; % Sampling period
T = 2000; % Final time (in seconds)
t = 0:dt:T-dt; % Time vector
L = length(t); % Length of the signal; # of samples

% Particle position data
x = zeros(1,L); % Particle's position for final SDE
x(1) = 0.5;   % Initial condition
y = x; % Particle's position for white noise

% Potential data
V0 = @(y) -y.^2./2 + y.^4./4; % Potential
dV0 = @(y) -y + y^3;
DV0 = V0(0) - V0(1); % Barrier height
V1 =  DV0*0.2; % Small compared to barrier
F = 1/100; % Frequency of the modulated potential
df1 = @(s) V1*sin(2*pi*F*s); % Modulated potential

kappa = .5; % Noise strength

M = 50; % Number of realisations
xM = zeros(M,L); % Mean of x
yM = xM; % Mean of y

for m = 1:M
    % EULER MARUYAMA method
    for n = 1:L-1
    % POTENTIAL NOISE: no té pics, però és semblant al POTENTIAL
        y(n+1) = y(n) - dt*dV0(y(n)) + sqrt(kappa*dt)*randn;
    % POTENTIAL PERIODIC NOISE
        x(n+1) = x(n) - dt*( dV0(x(n)) + df1(t(n+1)) ) + sqrt(kappa*dt)*randn; 
    end 
    xM(m,:) = x;
    yM(m,:) = y;
end
xmean = mean(xM,1);
ymean = mean(yM,1);
%%
subplot(4,2,1)
    plot(t,x)
    title('Single realisation of the final SDE','Interpreter','latex')
    xlabel('Time (s)')
    ylabel('Position x(t)')
subplot(4,2,2)
    periodogram(x,rectwin(L),L,Fs)
    xlim([0,.1])
subplot(4,2,3)
    plot(t,xmean)
    title(['Average of M = ',num2str(M),' realisations of the final SDE'],'Interpreter','latex')
    xlabel('Time (s)')
    ylabel('Position x(t)')
subplot(4,2,4)
    periodogram(xmean,rectwin(L),L,Fs)
    xlim([0,.1])
subplot(4,2,5)
    plot(t,y)
    title('Single realisation of the (only) noise SDE','Interpreter','latex')
    xlabel('Time (s)')
    ylabel('Position x(t)')
subplot(4,2,6)
    periodogram(y,rectwin(L),L,Fs)
    xlim([0,.1])
subplot(4,2,7)
    plot(t,ymean)
    title(['Average of M = ',num2str(M),' realisations of the (only) noise SDE'],'Interpreter','latex')
    xlabel('Time (s)')
    ylabel('Position x(t)')
subplot(4,2,8)
    periodogram(ymean,rectwin(L),L,Fs)
    xlim([0,.1])
sgtitle('Position of the particle and its periodogram','Interpreter', 'latex','FontSize',20)