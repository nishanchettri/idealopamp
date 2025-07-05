% Ideal_OpAmp_Generator
clear all;
close all;
clc;

% Define symbolic variables with real assumption
syms R1 R2 W P1 real

% Given constants
C1 = 1E-9;   % Farads
C2 = 1E-9;  % FaradsP2 (Hz)
P2_Hz = 20e+6; % Desired P2 frequency in Hz, try 2x the desired unity gain frequency
P2 = 2 * pi * P2_Hz; % Convert P2 to radians/sec
%P2 = 1.5e6;
A = 50; % DC Gain
PM = 60;    % Phase Margin in degrees

% Define numerator and denominator terms for symbolic equation
num = A * P1 * P2;
den = sqrt((W * P2)^2 + (10 * P2 - W^2)^2);

% Define equations to solve for P1 and W
eq1 = den == num;
eq2 = -atan((W * P2) / (P1 * P2 - W^2)) == PM * pi / 180;

% Solve equations for P1 and W with real constraints
S = solve([eq1, eq2], [P1, W], 'Real', true, 'IgnoreAnalyticConstraints', true);

P1_vals = double(S.P1);
W_vals = double(S.W);

% Filter for positive real solutions
valid_idx = find(imag(P1_vals) == 0 & imag(W_vals) == 0 & P1_vals > 0 & W_vals > 0);
if isempty(valid_idx)
    error('No valid real positive solutions found for P1 and W.');
end

P1 = P1_vals(valid_idx(1));
W = W_vals(valid_idx(1));

% Calculate Unity Gain Frequency (UGF) in Hz
Calculated_UGF = sqrt(2*(sqrt(4*A^2*P1^2*P2^2 + P1^4 - 2*P1^2*P2^2 + P2^4) - P1^2 - P2^2)) / (4 * pi);

% Solve for R1 and R2 based on pole equations
eqns = [ ...
    P1 == -(sqrt(C1^2*R1^2 + 2*C1*C2*R1*(R1 - R2) + C2^2*(R1^2 + 2*R1*R2 + R2^2)) - C1*R1 - C2*(R1 + R2)) / (2*C1*C2*R1*R2), ...
    P2 ==  (sqrt(C1^2*R1^2 + 2*C1*C2*R1*(R1 - R2) + C2^2*(R1^2 + 2*R1*R2 + R2^2)) + C1*R1 + C2*(R1 + R2)) / (2*C1*C2*R1*R2) ...
];

S = solve(eqns, [R1, R2], 'Real', true, 'IgnoreAnalyticConstraints', true);

R1_vals = double(S.R1);
R2_vals = double(S.R2);

% Filter for positive real resistor values
valid_idx = find(imag(R1_vals) == 0 & imag(R2_vals) == 0 & R1_vals > 0 & R2_vals > 0);
if isempty(valid_idx)
    error('No valid positive real solutions found for R1 and R2.');
end

R1 = R1_vals(valid_idx(1));
R2 = R2_vals(valid_idx(1));

% Recalculate poles P1 and P2 with solved R1 and R2 (rad/s)
P1_calc = -(sqrt(C1^2*R1^2 + 2*C1*C2*R1*(R1 - R2) + C2^2*(R1^2 + 2*R1*R2 + R2^2)) - C1*R1 - C2*(R1 + R2)) / (2*C1*C2*R1*R2);
P2_calc =  (sqrt(C1^2*R1^2 + 2*C1*C2*R1*(R1 - R2) + C2^2*(R1^2 + 2*R1*R2 + R2^2)) + C1*R1 + C2*(R1 + R2)) / (2*C1*C2*R1*R2);

% Define transfer function numerator and denominator (rad/s)
num_tf = [A * P1 * P2];
den_tf = [1, P1 + P2, P1 * P2];

% Create transfer function
sys = tf(num_tf, den_tf);

% Frequency vector for gain plot (1 Hz to 30 MHz for better coverage)
f = logspace(0, log10(30e6), 1000);
w = 2 * pi * f; % Angular frequency

% Calculate frequency response
[mag, phase] = bode(sys, w);
mag = squeeze(mag); % Magnitude (absolute)
phase = squeeze(phase); % Phase (degrees)

% Display calculated parameters
fprintf('Calculated Parameters:\n');
fprintf('----------------------\n');
fprintf('P1 (rad/s): %.4e\n', P1);
fprintf('P1 (Hz): %.4f\n', P1 / (2 * pi));
fprintf('P2 (rad/s): %.4e\n', P2);
fprintf('P2 (Hz): %.4f\n', P2 / (2 * pi));
fprintf('R1 (Ohms): %.4f\n', R1);
fprintf('R2 (Ohms): %.4f\n', R2);
fprintf('Calculated Unity Gain Frequency (UGF) (MHz): %.4f\n', Calculated_UGF/1000000);

% Plot Gain (Magnitude in dB)
figure;
semilogx(f, 20*log10(mag), 'LineWidth', 2);
grid on;
title('Gain Plot');
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
xlim([min(f) max(f)]);

% Plot Phase
figure;
semilogx(f, phase, 'LineWidth', 2);
grid on;
title('Phase Plot');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
xlim([min(f) max(f)]);

% Plot Bode plot using MATLAB built-in function
figure;
bode(sys);
grid on;
title('Bode Plot of the Ideal Op-Amp Transfer Function');
