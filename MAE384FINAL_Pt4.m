%%MAE384FinalProject
%%Final Project for class MAE384

%% Part 1 
clc;clear

% Time parameters
h = 1;    % Step size (1 day)
t0 = 0;   % Start of the simulation (day 0)
tf = 100; % End of the simulation (day 100)

% Beta and Gamma values
parameters = [2.0, 0.2;  % Measles
              1.0, 0.1;  % COVID-19
              0.3, 0.1]; % Influenza

% Population initial conditions
Total = 1000;   % Total population
S0 = 990;       % Initial susceptible population
I0 = 10;        % Initial infected population
R0 = 0;         % Initial recovered population

% Initialize a time vector
time = t0:h:tf; % Vector from 0 to 100 in step sizes of 1 day

% Initialize vector for plot titles
Disease = {'Measles', 'COVID-19', 'Influenza'};

% Loop for each disease
for k = 1:3
    beta = parameters (k, 1);
    gamma = parameters (k, 2);

    % Set up SIR state variables
    Sx = S0;
    Ix = I0;
    Rx = R0;

    % Initialize arrays and their first values
    Susceptible_1(1) = S0;
    Infected_1(1) = I0;
    Recovered_1(1) = R0;

    for i = 2:length(time)
    % calculate each of the k values for for the SIR model
    [Sk1, Ik1, Rk1] = dSIRdt (Sx, Ix, Rx, beta, gamma, Total); 
    [Sk2, Ik2, Rk2] = dSIRdt (Sx + (.5 * Sk1 * h), ... % k2 susceptible
                              Ix + (.5 * Ik1 * h), ... % k2 infected
                              Rx + (.5 * Rk1 * h), ... % k2 recovered
                              beta, gamma, Total); % respective beta and gamma values
    [Sk3, Ik3, Rk3] = dSIRdt (Sx + (.5 * Sk2 * h), ... % k3 susceptible
                              Ix + (.5 * Ik2 * h), ... % k3 infected
                              Rx + (.5 * Rk2 * h), ... % k3 recovered
                              beta, gamma, Total); % respective beta and gamma values
    [Sk4, Ik4, Rk4] = dSIRdt (Sx + (Sk3 * h), ... % k4 susceptible
                              Ix + (Ik3 * h), ... % k4 infected
                              Rx + (Rk3 * h), ... % k4 recovered
                              beta, gamma, Total); % respective beta and gamma values

    % calculating the new values via the Runge-Kutta formula
    Sx = Sx + ((h / 6) * (Sk1 + 2*Sk2 + 2*Sk3 + Sk4)); % RK susceptible
    Ix = Ix + ((h / 6) * (Ik1 + 2*Ik2 + 2*Ik3 + Ik4)); % RK infected
    Rx = Rx + ((h / 6) * (Rk1 + 2*Rk2 + 2*Rk3 + Rk4)); % RK recovered

    % Store calculated values
    Susceptible_1(i) = Sx;
    Infected_1(i) = Ix;
    Recovered_1(i) = Rx;

    end

% Create a plot for each of the groups

figure;
plot (time, Susceptible_1,'r', time, Infected_1,'g', time, Recovered_1,'b')
legend('Susceptible', 'Infected', 'Recovered')
xlabel ('Time (days)')
ylabel ('Respective Populations')
title (Disease{k})
grid on

end



%% Part 2
% Calculate SIR model for step size of 2 days for influenza
% Time parameters
h2 = 2;    % Step size (2 days)
t0 = 0;   % Start of the simulation (day 0)
tf = 100; % End of the simulation (day 100)

% Initialize a time vector
time = t0:h2:tf; % Vector from 0 to 100 in step sizes of 1 day

    % Set up SIR state variables
    Sx = S0;
    Ix = I0;
    Rx = R0;

    % Initialize arrays and their first values
    Susceptible_even(1) = S0;
    Infected_even(1) = I0;
    Recovered_even(1) = R0;

    for i = 2:length(time)
    % calculate each of the k values for for the SIR model
    [Sk1, Ik1, Rk1] = dSIRdt (Sx, Ix, Rx, beta, gamma, Total); 
    [Sk2, Ik2, Rk2] = dSIRdt (Sx + (.5 * Sk1 * h2), ... % k2 susceptible
                              Ix + (.5 * Ik1 * h2), ... % k2 infected
                              Rx + (.5 * Rk1 * h2), ... % k2 recovered
                              beta, gamma, Total); % respective beta and gamma values
    [Sk3, Ik3, Rk3] = dSIRdt (Sx + (.5 * Sk2 * h2), ... % k3 susceptible
                              Ix + (.5 * Ik2 * h2), ... % k3 infected
                              Rx + (.5 * Rk2 * h2), ... % k3 recovered
                              beta, gamma, Total); % respective beta and gamma values
    [Sk4, Ik4, Rk4] = dSIRdt (Sx + (Sk3 * h2), ... % k4 susceptible
                              Ix + (Ik3 * h2), ... % k4 infected
                              Rx + (Rk3 * h2), ... % k4 recovered
                              beta, gamma, Total); % respective beta and gamma values

    % calculating the new values via the Runge-Kutta formula
    Sx = Sx + ((h2 / 6) * (Sk1 + 2*Sk2 + 2*Sk3 + Sk4)); % RK susceptible
    Ix = Ix + ((h2 / 6) * (Ik1 + 2*Ik2 + 2*Ik3 + Ik4)); % RK infected
    Rx = Rx + ((h2 / 6) * (Rk1 + 2*Rk2 + 2*Rk3 + Rk4)); % RK recovered

    % Store calculated values
    Susceptible_even(i) = Sx;
    Infected_even(i) = Ix;
    Recovered_even(i) = Rx;

    end

%% Set up for interpolation

% Define given data (e.g., from the SIR model with coarse time step h = 2)
t_even = 0:2:100;       % Time points where data is available (even days)
S_even = Susceptible_even; % Susceptible values at t_even
I_even = Infected_even;    % Infected values at t_even
R_even = Recovered_even;   % Recovered values at t_even

% Define the real values with the smaller stepsizes from part 1
S_real = Susceptible_1(1:2:99); % Susceptible values at t odd, real
I_real = Infected_1(1:2:99);    % Infected values at t odd, real
R_real = Recovered_1(1:2:99);   % Recovered values at t odd, real

% Define time points for interpolation (odd days)
t_odd = 1:2:99;

% Initialize interpolated values for linear and quadratic
S_interp = zeros(size(t_odd));
I_interp = zeros(size(t_odd));
R_interp = zeros(size(t_odd));

Linear_S_interp = zeros(size(t_odd));
Linear_I_interp = zeros(size(t_odd));
Linear_R_interp = zeros(size(t_odd));

%% Linear interpolation using Lagrange functions

lagrange_linear = @(x, x1, y1, x2, y2) ...
    y1 * ((x - x2) / (x1 - x2)) + y2 * ((x - x1) / (x2 - x1));

for k = 1:length(t_odd)
    t_k = t_odd(k);

    % Find the first point in t_even greater than t_k
    idx = find(t_even > t_k, 1); 

    % Handle boundaries: ensure idx-1 and idx+1 are valid
    if idx == 1
        % Use the first two points if idx is at the start
        Lx1 = t_even(idx); 
        Lx2 = t_even(idx+1); 

        Ly1_S = S_even(idx); 
        Ly2_S = S_even(idx+1); 

        Ly1_I = I_even(idx); 
        Ly2_I = I_even(idx+1); 

        Ly1_R = R_even(idx); 
        Ly2_R = R_even(idx+1); 
    elseif idx == length(t_even)
        % Use the last three points if idx is at the end
        Lx1 = t_even(idx-1); 
        Lx2 = t_even(idx);

        Ly1_S = S_even(idx-1); 
        Ly2_S = S_even(idx);

        Ly1_I = I_even(idx-1); 
        Ly2_I = I_even(idx);

        Ly1_R = R_even(idx-1); 
        Ly2_R = R_even(idx);
    else
        % Normal case: take forward points 
        Lx1 = t_even(idx); 
        Lx2 = t_even(idx+1);

        Ly1_S = S_even(idx); 
        Ly2_S = S_even(idx+1);

        Ly1_I = I_even(idx); 
        Ly2_I = I_even(idx+1);

        Ly1_R = R_even(idx); 
        Ly2_R = R_even(idx+1);
    end


    % Interpolate using the Lagrange formula
    Linear_S_interp(k) = lagrange_linear(t_k, Lx1, Ly1_S, Lx2, Ly2_S);
    Linear_I_interp(k) = lagrange_linear(t_k, Lx1, Ly1_I, Lx2, Ly2_I);
    Linear_R_interp(k) = lagrange_linear(t_k, Lx1, Ly1_R, Lx2, Ly2_R);
end

even = 1;
odd = 1;
for i = 1:length(Susceptible_1)
    if (-1)^i < 0
        S_linear(i) = S_even(even);
        I_linear(i) = I_even(even);
        R_linear(i) = R_even(even);
        even = even + 1;
    end 
    if (-1)^i > 0
        S_linear(i) = Linear_S_interp(odd);
        I_linear(i) = Linear_I_interp(odd);
        R_linear(i) = Linear_R_interp(odd);
        odd = odd + 1;
    end
end

Linear_S_error = sqrt((sum(S_linear-Susceptible_1).^2)./50);
Linear_I_error = sqrt((sum(I_linear-Infected_1).^2)./50);
Linear_R_error = sqrt((sum(R_linear-Recovered_1).^2)./50);


%% Quadratic interpolation using Lagrange functions

% Quadratic Lagrange interpolation function
lagrange_quadratic = @(x, x1, y1, x2, y2, x3, y3) ...
    y1 * ((x - x2) * (x - x3)) / ((x1 - x2) * (x1 - x3)) + ...
    y2 * ((x - x1) * (x - x3)) / ((x2 - x1) * (x2 - x3)) + ...
    y3 * ((x - x1) * (x - x2)) / ((x3 - x1) * (x3 - x2));

for k = 1:length(t_odd)
    t_k = t_odd(k);

    % Find the first point in t_even greater than t_k
    idx = find(t_even > t_k, 1); 

    % Handle boundaries: ensure idx-1 and idx+1 are valid
    if idx == 1
        % Use the first three points if idx is at the start
        x1 = t_even(idx); 
        x2 = t_even(idx+1); 
        x3 = t_even(idx+2);

        y1_S = S_even(idx); 
        y2_S = S_even(idx+1); 
        y3_S = S_even(idx+2);

        y1_I = I_even(idx); 
        y2_I = I_even(idx+1); 
        y3_I = I_even(idx+2);

        y1_R = R_even(idx); 
        y2_R = R_even(idx+1); 
        y3_R = R_even(idx+2);
    elseif idx == length(t_even)
        % Use the last three points if idx is at the end
        x1 = t_even(idx-2); 
        x2 = t_even(idx-1); 
        x3 = t_even(idx);

        y1_S = S_even(idx-2); 
        y2_S = S_even(idx-1); 
        y3_S = S_even(idx);

        y1_I = I_even(idx-2); 
        y2_I = I_even(idx-1); 
        y3_I = I_even(idx);

        y1_R = R_even(idx-2); 
        y2_R = R_even(idx-1); 
        y3_R = R_even(idx);
    else
        % Normal case: take three nearest points
        x1 = t_even(idx-1); 
        x2 = t_even(idx); 
        x3 = t_even(idx+1);

        y1_S = S_even(idx-1); 
        y2_S = S_even(idx); 
        y3_S = S_even(idx+1);

        y1_I = I_even(idx-1); 
        y2_I = I_even(idx); 
        y3_I = I_even(idx+1);

        y1_R = R_even(idx-1); 
        y2_R = R_even(idx); 
        y3_R = R_even(idx+1);
    end

    % Interpolate using the Lagrange formula
    S_interp(k) = lagrange_quadratic(t_k, x1, y1_S, x2, y2_S, x3, y3_S);
    I_interp(k) = lagrange_quadratic(t_k, x1, y1_I, x2, y2_I, x3, y3_I);
    R_interp(k) = lagrange_quadratic(t_k, x1, y1_R, x2, y2_R, x3, y3_R);

end

even = 1;
odd = 1;
for i = 1:length(Susceptible_1)
    if (-1)^i < 0
        S_quad(i) = S_even(even);
        I_quad(i) = I_even(even);
        R_quad(i) = R_even(even);
        even = even + 1;
    end 
    if (-1)^i > 0
        S_quad(i) = S_interp(odd);
        I_quad(i) = I_interp(odd);
        R_quad(i) = R_interp(odd);
        odd = odd + 1;
    end
end
% Compute error in each part
S_error = sqrt((sum(S_quad-Susceptible_1).^2)./50);
I_error = sqrt((sum(I_quad-Infected_1).^2)./50);
R_error = sqrt((sum(R_quad-Recovered_1).^2)./50);

Error = [Linear_S_error, S_error;
         Linear_I_error, I_error; 
         Linear_R_error, R_error];

% Display Error
disp('Error (Linear, Quadratic):');
disp(Error)
disp('Rows are S, I, R, respectively')

% Discussion: Quadratic Lagrange interpolation results in smaller errors.
% This makes sense because it uses more points and thus is more accurate.

%% Part 3: Least Squares Estimation

% Run the non-linear SIR model from Part 1 to generate "true" I(t)
h = 1; % Step size
t0 = 0; % Start of simulation
tf = 30; % End of simulation (30 days)
time = t0:h:tf; % Time vector
Total = 1000; % Total population
S0 = 990; % Initial susceptible population
I0_true = 10; % Initial infected population
R0 = 0; % Initial recovered population
beta_true = 0.3; % Transmission rate
gamma_true = 0.1; % Recovery rate

% Initialize state variables
S = S0;
I = I0_true;
R = R0;
I_true = zeros(size(time));
I_true(1) = I0_true;

% Run the non-linear SIR model
for i = 2:length(time)
    [dS, dI, dR] = dSIRdt(S, I, R, beta_true, gamma_true, Total);
    S = S + h * dS;
    I = I + h * dI;
    R = R + h * dR;
    I_true(i) = I;
end

% Take the natural logarithm of I(t) to linearize
ln_I = log(I_true);

% Apply linear least squares to estimate k and I(0)
X = [ones(length(time), 1), time']; % Design matrix: [1, t]
Y = ln_I'; % Observed data
coeff = (X \ Y); % Least squares solution

ln_I0_est = coeff(1); % Estimate of ln(I(0))
k_est = coeff(2); % Estimate of k

% Convert ln(I0) to I0
I0_est = exp(ln_I0_est);

% Estimate beta from k
k_true = (beta_true * S0 / Total) - gamma_true; % True k value
beta_est = (k_est + gamma_true) * Total / S0;

fprintf('Results for 30 Days of Data:\n');
fprintf('Estimated I(0): %.4f (True I(0): %.4f)\n', I0_est, I0_true);
fprintf('Estimated k: %.4f (True k: %.4f)\n', k_est, k_true);
fprintf('Estimated beta: %.4f (True beta: %.4f)\n\n', beta_est, beta_true);

%% Repeat with 10 Days of Data
time_10 = time(1:10); % Use first 10 days of data
ln_I_10 = ln_I(1:10); % Corresponding ln(I(t))

% Apply linear least squares to the reduced dataset
X_10 = [ones(length(time_10), 1), time_10'];
Y_10 = ln_I_10';
coeff_10 = (X_10 \ Y_10);

ln_I0_est_10 = coeff_10(1);
k_est_10 = coeff_10(2);

% Convert ln(I0) to I0 for 10 days
I0_est_10 = exp(ln_I0_est_10);

% Estimate beta from k for 10 days
beta_est_10 = (k_est_10 + gamma_true) * Total / S0;

fprintf('Results for 10 Days of Data:\n');
fprintf('Estimated I(0): %.4f (True I(0): %.4f)\n', I0_est_10, I0_true);
fprintf('Estimated k: %.4f (True k: %.4f)\n', k_est_10, k_true);
fprintf('Estimated beta: %.4f (True beta: %.4f)\n\n', beta_est_10, beta_true);

% Analysis
fprintf('Analysis:\n');
fprintf('- Using 30 days of data results in more accurate estimates compared to 10 days.\n');
fprintf('- With fewer data points (10 days), the error increases due to less information.\n');
fprintf('- Since the Part-1 graphs show a parabola the quadradic method is more accurate as opposed to a linear method.\n');

%% Part 4 : Fourier Analysis

%%Part 1 omega = 2*pi*(365/365)

% Time parameters
h = 0.1;    % Step size (0.1 day)
t0 = 0;   % Start of the simulation (day 0)
tf = 30; % End of the simulation (day 30)

% Initialize a time vector
time = t0:h:tf; % Vector from 0 to 30 in step sizes of 0.1 day

% Beta and Gamma values
Beta = 0.3.*(1+5.*sin(2.*pi.*(365/365).*time));
Gamma = 0.1;

% Population initial conditions
Total = 1000;   % Total population
S0 = 990;       % Initial susceptible population
I0 = 10;        % Initial infected population
R0 = 0;         % Initial recovered population

% Set up SIR state variables
Sx = S0;
Ix = I0;
Rx = R0;

% Initialize arrays and their first values
Susceptible_1 = zeros(1,length(time));
Infected_1 = zeros(1,length(time));
Recovered_1 = zeros(1,length(time));
Susceptible_1(1) = S0;
Infected_1(1) = I0;
Recovered_1(1) = R0;

for i = 2:length(time)
    % calculate each of the k values for for the SIR model
    [Sk1, Ik1, Rk1] = dSIRdt (Sx, Ix, Rx, Beta(i), Gamma, Total); 
    [Sk2, Ik2, Rk2] = dSIRdt (Sx + (.5 * Sk1 * h), ... % k2 susceptible
                              Ix + (.5 * Ik1 * h), ... % k2 infected
                              Rx + (.5 * Rk1 * h), ... % k2 recovered
                              Beta(i), Gamma, Total); % respective beta and gamma values
    [Sk3, Ik3, Rk3] = dSIRdt (Sx + (.5 * Sk2 * h), ... % k3 susceptible
                              Ix + (.5 * Ik2 * h), ... % k3 infected
                              Rx + (.5 * Rk2 * h), ... % k3 recovered
                              Beta(i), Gamma, Total); % respective beta and gamma values
    [Sk4, Ik4, Rk4] = dSIRdt (Sx + (Sk3 * h), ... % k4 susceptible
                              Ix + (Ik3 * h), ... % k4 infected
                              Rx + (Rk3 * h), ... % k4 recovered
                              Beta(i), Gamma, Total); % respective beta and gamma values
    
    % calculating the new values via the Runge-Kutta formula
    Sx = Sx + ((h / 6) * (Sk1 + 2*Sk2 + 2*Sk3 + Sk4)); % RK susceptible
    Ix = Ix + ((h / 6) * (Ik1 + 2*Ik2 + 2*Ik3 + Ik4)); % RK infected
    Rx = Rx + ((h / 6) * (Rk1 + 2*Rk2 + 2*Rk3 + Rk4)); % RK recovered
    
    % Store calculated values
    Susceptible_1(i) = Sx;
    Infected_1(i) = Ix;
    Recovered_1(i) = Rx;

end

% Create a plot for first Beta functiuon

figure;
plot (time, Susceptible_1,'r', time, Infected_1,'g', time, Recovered_1,'b')
legend('Susceptible', 'Infected', 'Recovered')
xlabel ('Time')
ylabel ('Respective Populations')
title ('Variable Spread Rate Omega = 2*Pi*(365/365)')
grid on


% Define frequency range
T = 30; % Total length of signal (30 days)
N = length(time); % Number of samples
f = (1/T)*(0:(N/2)); % Frequency Vector

% Calculate spectrum
spectrum = fft(Infected_1);

P2 = abs(spectrum/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure;
plot(f,P1);
title("Spectrum of Infected Cases, Omega = 2*Pi*(365/365)")
xlabel("Frequency Coefficients")
ylabel("|fft(I(t))|")



%%Part 2 omega = 2*pi*(100/365)

clearvars;


% Time parameters
h = 0.1;    % Step size (0.1 day)
t0 = 0;   % Start of the simulation (day 0)
tf = 30; % End of the simulation (day 30)

% Initialize a time vector
time = t0:h:tf; % Vector from 0 to 30 in step sizes of 0.1 day

% Beta and Gamma values
Beta = 0.3.*(1+5.*sin(2.*pi.*(100/365).*time));
Gamma = 0.1;

% Population initial conditions
Total = 1000;   % Total population
S0 = 990;       % Initial susceptible population
I0 = 10;        % Initial infected population
R0 = 0;         % Initial recovered population

% Set up SIR state variables
Sx = S0;
Ix = I0;
Rx = R0;

% Initialize arrays and their first values
Susceptible_1 = zeros(1,length(time));
Infected_1 = zeros(1,length(time));
Recovered_1 = zeros(1,length(time));
Susceptible_1(1) = S0;
Infected_1(1) = I0;
Recovered_1(1) = R0;

for i = 2:length(time)
    % calculate each of the k values for for the SIR model
    [Sk1, Ik1, Rk1] = dSIRdt (Sx, Ix, Rx, Beta(i), Gamma, Total); 
    [Sk2, Ik2, Rk2] = dSIRdt (Sx + (.5 * Sk1 * h), ... % k2 susceptible
                              Ix + (.5 * Ik1 * h), ... % k2 infected
                              Rx + (.5 * Rk1 * h), ... % k2 recovered
                              Beta(i), Gamma, Total); % respective beta and gamma values
    [Sk3, Ik3, Rk3] = dSIRdt (Sx + (.5 * Sk2 * h), ... % k3 susceptible
                              Ix + (.5 * Ik2 * h), ... % k3 infected
                              Rx + (.5 * Rk2 * h), ... % k3 recovered
                              Beta(i), Gamma, Total); % respective beta and gamma values
    [Sk4, Ik4, Rk4] = dSIRdt (Sx + (Sk3 * h), ... % k4 susceptible
                              Ix + (Ik3 * h), ... % k4 infected
                              Rx + (Rk3 * h), ... % k4 recovered
                              Beta(i), Gamma, Total); % respective beta and gamma values
    
    % calculating the new values via the Runge-Kutta formula
    Sx = Sx + ((h / 6) * (Sk1 + 2*Sk2 + 2*Sk3 + Sk4)); % RK susceptible
    Ix = Ix + ((h / 6) * (Ik1 + 2*Ik2 + 2*Ik3 + Ik4)); % RK infected
    Rx = Rx + ((h / 6) * (Rk1 + 2*Rk2 + 2*Rk3 + Rk4)); % RK recovered
    
    % Store calculated values
    Susceptible_1(i) = Sx;
    Infected_1(i) = Ix;
    Recovered_1(i) = Rx;

end

% Create a plot for first Beta functiuon

figure;
plot (time, Susceptible_1,'r', time, Infected_1,'g', time, Recovered_1,'b')
legend('Susceptible', 'Infected', 'Recovered')
xlabel ('Time')
ylabel ('Respective Populations')
title ('Variable Spread Rate Omega = 2*Pi*(100/365)')
grid on


% Define frequency range
T = 30; % Total length of signal (30 days)
N = length(time); % Number of samples
f = (1/T)*(0:(N/2)); % Frequency Vector

% Calculate spectrum
spectrum = fft(Infected_1);

P2 = abs(spectrum/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure;
plot(f,P1);
title("Spectrum of Infected Cases, Omega = 2*Pi*(100/365)")
xlabel("Frequency Coefficients")
ylabel("|fft(I(t))|")

% Define function for derivatives (Or nothing will work)

function [dS, dI, dR] = dSIRdt(S, I, R, beta, gamma, N)

    % Create a function that calculates the derivatives of the susceptible,
    % infected, and recovered population with respect of time (dS/dt,
    % dI/dt, dR/dt).
    % S = Susceptible population
    % I = Infected population
    % R = Recovered population
    % beta = transmission rate (how often a susceptible individual gets
    % infected)
    % gamma = recovery rate (how often an infected individual recovers)
    % N = total population (S + I + R = constant)

    dS = -(beta / N) * S * I; % define the derivative of the susceptible population
    dI = ((beta / N) * S * I) - (gamma * I); % Define the derivative of the infected population
    dR = gamma * I; % Define the derivative of the infected population
end