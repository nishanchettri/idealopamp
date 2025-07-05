% Ideal_OpAmp_Generator
clear all;
close all;
clc;

syms R1 R2 W P1

C1 = 1E-6;
C2 = 1E-6;

P2 = 1.5E6; % Radians/sec
A = 100000; % DC Gain
PM = 60;    % Phase Margin in degrees

% Solve for P1 and W from equations
num = [A*P1*P2];
den = sqrt((W*P2)^2 + (10*P2 - W^2)^2);

eqs = [den == num, -atan((W*P2/(P1*P2 - W^2))) == PM*pi/180];
S = solve(eqs, [P1 W]);
P1 = double(S.P1);
W = double(S.W);

% Convert from radians/sec to Hz
P1 = P1(1)/(2*pi);
P2 = P2/(2*pi);

% Calculate Unity Gain Frequency (UGF)
Calculated_UGF = sqrt(2*(sqrt(4*A^2*P1^2*P2^2 + P1^4 - 2*P1^2*P2^2 + P2^4) - P1^2 - P2^2))/(4*pi);

% Solve for R1 and R2
eqns = [P1 == -(sqrt((C1)^2*(R1)^2 + 2*C1*C2*R1*(R1-R2) + (C2)^2*((R1)^2 + 2*R1*R2 + (R2)^2)) - C1*R1 - C2*(R1 + R2))/(2*C1*C2*R1*R2), 
        P2 == (sqrt((C1)^2*(R1)^2 + 2*C1*C2*R1*(R1-R2) + (C2)^2*((R1)^2 + 2*R1*R2 + (R2)^2)) + C1*R1 + C2*(R1 + R2))/(2*C1*C2*R1*R2)];
S = solve(eqns, [R1 R2]);
R1 = double(S.R1(1));
R2 = double(S.R2(1));

% Recalculate poles P1 and P2 with solved R1 and R2
P1 = -(sqrt((C1)^2*(R1)^2 + 2*C1*C2*R1*(R1-R2) + (C2)^2*((R1)^2 + 2*R1*R2 + (R2)^2)) - C1*R1 - C2*(R1 + R2)) / (2*C1*C2*R1*R2);
P2 = (sqrt((C1)^2*(R1)^2 + 2*C1*C2*R1*(R1-R2) + (C2)^2*((R1)^2 + 2*R1*R2 + (R2)^2)) + C1*R1 + C2*(R1 + R2)) / (2*C1*C2*R1*R2);

% Define transfer function numerator and denominator
num = [A*P1*P2];
denom = [1, P1 + P2, P1*P2];

% Create transfer function
sys = tf(num, denom);

% Frequency vector for gain plot
f = logspace(0, 7, 1000); % Frequency from 1 Hz to 10 MHz
w = 2*pi*f; % Angular frequency

% Calculate frequency response
[mag, phase] = bode(sys, w);
mag = squeeze(mag); % Magnitude in absolute
phase = squeeze(phase); % Phase in degrees

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

% After all calculations, add:

fprintf('Calculated Parameters:\n');
fprintf('----------------------\n');
fprintf('P1 (rad/s): %.4e\n', P1);
fprintf('P1 (Hz): %.4f\n', P1/(2*pi));
fprintf('P2 (rad/s): %.4e\n', P2*2*pi); % Remember you converted P2 to Hz earlier
fprintf('P2 (Hz): %.4f\n', P2);
fprintf('W (rad/s): %.4e\n', W);
fprintf('W (Hz): %.4f\n', W/(2*pi));
fprintf('R1 (Ohms): %.4f\n', R1);
fprintf('R2 (Ohms): %.4f\n', R2);
fprintf('Calculated Unity Gain Frequency (UGF) (Hz): %.4f\n', Calculated_UGF);
