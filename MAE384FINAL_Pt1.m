% MAE384FinalProject
% Final Project for class MAE384

clc;clear

% Time parameters
h = 1;    % Step size (1 day)
t0 = 0;   % Start of the simulation (day 0)
tf = 100; % End of the simulation (day 100)

% Beta and Gamma values
parameters = [0.3, 0.1;  % Influenza
              1.0, 0.1;  % COVID-19
              2.0, 0.2]; % Measles

% Population initial conditions
Total = 1000;   % Total population
S0 = 990;       % Initial susceptible population
I0 = 10;        % Initial infected population
R0 = 0;         % Initial recovered population

% Initialize a time vector
time = t0:h:tf; % Vector from 0 to 100 in step sizes of 1 day

% Initialize vector for plot titles
Disease = {'Influenza', 'COVID-19', 'Measles'};

% Loop for each disease
for k = 1:3
    beta = parameters (k, 1);
    gamma = parameters (k, 2);
    
    % Set up SIR state variables
    Sx = S0;
    Ix = I0;
    Rx = R0;

    % Initialize arrays and their first values
    Susceptible(1) = S0;
    Infected(1) = I0;
    Recovered(1) = R0;

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
    Susceptible(i) = Sx;
    Infected(i) = Ix;
    Recovered(i) = Rx;

    end

% Create a plot for each of the groups

figure;
plot (time, Susceptible,'r', time, Infected,'g', time, Recovered,'b')
legend('Susceptible', 'Infected', 'Recovered')
xlabel ('Time (days)')
ylabel ('Respective Populations')
title (Disease{k})
grid on

end

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
