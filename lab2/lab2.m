% ELE 532 Signals & Systems I - Laboratory Assignment 2
% Final Script - English Only

% --- Environment Initialization ---
clear all;
close all;
clc;

%% ========================================================================
%  Create Main Figure and Tab Group
% =========================================================================

mainFig = figure('Name', 'ELE 532 Lab 2 Results (All Problems)', 'NumberTitle', 'off', 'WindowState', 'maximized');
tabGroup = uitabgroup('Parent', mainFig);

%% ========================================================================
%  A. Impulse Response
% =========================================================================

fprintf('Executing Problem A.1 & A.2...\n');

% --- Problem A.1 & A.2: Characteristic Polynomial and Impulse Response Plot ---

lambda_A = [-5, -10];
p_A = poly(lambda_A);

fprintf('Problem A.1: The coefficients of the characteristic polynomial are:\n');
disp(p_A);

t_A = 0:0.0005:0.1;
h_A = exp(lambda_A(1)*t_A) - exp(lambda_A(2)*t_A);

tab_A2 = uitab('Parent', tabGroup, 'Title', 'Problem A.2');
ax_A2 = axes('Parent', tab_A2);

plot(ax_A2, t_A, h_A, 'LineWidth', 2, 'Color', 'b');
title(ax_A2, 'Problem A.2: System Impulse Response h(t)');
xlabel(ax_A2, 'Time (t)');
ylabel(ax_A2, 'Amplitude');
grid(ax_A2, 'on');
fprintf('Graph for Problem A.2 has been generated in its tab.\n\n');

%% ========================================================================
%  B. Convolution
% =========================================================================

% --- Problem B.2: Convolution of a Rectangular and a Triangular Pulse ---
fprintf('Executing Problem B.2...\n');

dt_B2 = 0.01;
t_B2 = -5:dt_B2:10;

x_B2 = (t_B2 >= 0 & t_B2 <= 3);

h_B2 = zeros(size(t_B2));
h_B2(t_B2 >= 0 & t_B2 <= 2) = t_B2(t_B2 >= 0 & t_B2 <= 2);
h_B2(t_B2 > 2 & t_B2 <= 4) = 4 - t_B2(t_B2 > 2 & t_B2 <= 4);

y_B2 = conv(x_B2, h_B2) * dt_B2;
t_start_y_B2 = t_B2(1) + t_B2(1);
t_y_B2 = t_start_y_B2:dt_B2:(t_start_y_B2 + (length(y_B2)-1)*dt_B2);

tab_B2 = uitab('Parent', tabGroup, 'Title', 'Problem B.2');

subplot(3,1,1, 'Parent', tab_B2); plot(t_B2, x_B2, 'LineWidth', 2); title('Input Signal x(t)'); grid on; xlabel('Time (t)'); ylabel('Amplitude');
subplot(3,1,2, 'Parent', tab_B2); plot(t_B2, h_B2, 'LineWidth', 2); title('Impulse Response h(t)'); grid on; xlabel('Time (t)'); ylabel('Amplitude');
subplot(3,1,3, 'Parent', tab_B2); plot(t_y_B2, y_B2, 'LineWidth', 2); title('Convolution Output y(t) = x(t) * h(t)'); grid on; xlabel('Time (t)'); ylabel('Amplitude');
fprintf('Graph for Problem B.2 has been generated in its tab.\n\n');

% --- Problem B.3: Convolution of two Rectangular Pulses ---
fprintf('Executing Problem B.3...\n');

dt_B3 = 0.01;
t_B3 = -2:dt_B3:8;

x1_B3 = (t_B3 >= 0 & t_B3 <= 4);
x2_B3 = (t_B3 >= 0 & t_B3 <= 2);

y_B3 = conv(x1_B3, x2_B3) * dt_B3;
t_start_y_B3 = t_B3(1) + t_B3(1);
t_y_B3 = t_start_y_B3:dt_B3:(t_start_y_B3 + (length(y_B3)-1)*dt_B3);

tab_B3 = uitab('Parent', tabGroup, 'Title', 'Problem B.3');

subplot(3,1,1, 'Parent', tab_B3); plot(t_B3, x1_B3, 'LineWidth', 2); title('Signal x1(t)'); grid on; ylim([0 1.2]);
subplot(3,1,2, 'Parent', tab_B3); plot(t_B3, x2_B3, 'LineWidth', 2); title('Signal x2(t)'); grid on; ylim([0 1.2]);
subplot(3,1,3, 'Parent', tab_B3); plot(t_y_B3, y_B3, 'LineWidth', 2); title('Convolution Output y(t) = x1(t) * x2(t)'); grid on;
fprintf('Graph for Problem B.3 has been generated in its tab.\n\n');

%% ========================================================================
%  C. System Behavior and Stability
% =========================================================================

% --- Problem C.1 & C.2: Impulse Response Plots and Eigenvalues ---
fprintf('Executing Problem C.1 & C.2...\n');

t_C1 = -1:0.001:5;
u_C1 = (t_C1 >= 0);

h1_C1 = exp(0.5 * t_C1) .* u_C1;
h2_C1 = 4 * exp(-0.5 * t_C1) .* u_C1;
h3_C1 = 4 * exp(-t_C1) .* u_C1;
h4_C1 = 4 * (exp(-0.5 * t_C1) - exp(-t_C1)) .* u_C1;

tab_C1 = uitab('Parent', tabGroup, 'Title', 'Problem C.1');
ax_C1 = axes('Parent', tab_C1);

plot(ax_C1, t_C1, h1_C1, 'LineWidth', 2); hold(ax_C1, 'on');
plot(ax_C1, t_C1, h2_C1, 'LineWidth', 2);
plot(ax_C1, t_C1, h3_C1, 'LineWidth', 2);
plot(ax_C1, t_C1, h4_C1, 'LineWidth', 2); hold(ax_C1, 'off');
title(ax_C1, 'Problem C.1: Impulse Responses of Systems S1-S4');
xlabel(ax_C1, 'Time (t)');
ylabel(ax_C1, 'Amplitude');
legend(ax_C1, 'h1(t) - Unstable', 'h2(t) - Stable', 'h3(t) - Stable', 'h4(t) - Stable');
grid(ax_C1, 'on');
axis(ax_C1, [-1 5 -1 5]);
fprintf('Graph for Problem C.1 has been generated in its tab.\n');

fprintf('\nProblem C.2: Characteristic values (Eigenvalues) of the systems:\n');
fprintf('  - S1: %.1f (Real part is positive => Unstable)\n', 0.5);
fprintf('  - S2: %.1f (Real part is negative => Stable)\n', -0.5);
fprintf('  - S3: %.1f (Real part is negative => Stable)\n', -1.0);
fprintf('  - S4: %.1f, %.1f (Real parts are negative => Stable)\n\n', -0.5, -1.0);

% --- Problem C.3: Convolution with Truncated Impulse Responses ---
fprintf('Executing Problem C.3...\n');

dtau_C3 = 0.01;
tau_C3 = 0:dtau_C3:20;
u_tau_C3 = (tau_C3 >= 0);

h1_trunc = exp(0.5 * tau_C3) .* u_tau_C3;
h2_trunc = 4 * exp(-0.5 * tau_C3) .* u_tau_C3;
h3_trunc = 4 * exp(-tau_C3) .* u_tau_C3;
h4_trunc = 4 * (exp(-0.5 * tau_C3) - exp(-tau_C3)) .* u_tau_C3;

tvec_C3 = 0:0.1:20;
x_input_C3 = (sin(5 * tvec_C3)) .* (tvec_C3 >= 0 & tvec_C3 < 3);

y1_out = conv(h1_trunc, x_input_C3) * dtau_C3;
y2_out = conv(h2_trunc, x_input_C3) * dtau_C3;
y3_out = conv(h3_trunc, x_input_C3) * dtau_C3;
y4_out = conv(h4_trunc, x_input_C3) * dtau_C3;

t_start_out = tau_C3(1) + tvec_C3(1);
t_out_length = length(y1_out);
t_out = t_start_out:dtau_C3:(t_start_out + (t_out_length-1)*dtau_C3);

tab_C3 = uitab('Parent', tabGroup, 'Title', 'Problem C.3');
ax_C3 = axes('Parent', tab_C3);

plot(ax_C3, t_out, y1_out, 'DisplayName', 'Output of S1 (Unstable)'); hold(ax_C3, 'on');
plot(ax_C3, t_out, y2_out, 'DisplayName', 'Output of S2');
plot(ax_C3, t_out, y3_out, 'DisplayName', 'Output of S3');
plot(ax_C3, t_out, y4_out, 'DisplayName', 'Output of S4'); hold(ax_C3, 'off');
title(ax_C3, 'Problem C.3: System Outputs for a Given Input');
xlabel(ax_C3, 'Time (t)');
ylabel(ax_C3, 'Output Amplitude');
legend(ax_C3, 'show');
grid(ax_C3, 'on');
fprintf('Graph for Problem C.3 has been generated in its tab.\n');

fprintf('\nAll problems have been executed. Check the figure window with tabs.\n');