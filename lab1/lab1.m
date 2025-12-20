% Name: keunhyeok choi
% Student ID: 500958125
% ELE 532: Signals & Systems I
% Laboratory Assignment 1

clear;
clc;
close all;

%% =====================================================================
% A. Anonymous Functions and Plotting
% ======================================================================
disp('Executing Section A: Anonymous Functions and Plotting...');

% Problem A.1
f_exp = @(t) exp(-t);
t1 = -2:0.01:2;
figure('Name','Section A Plots');
subplot(2, 1, 1);
plot(t1, f_exp(t1));
title('Problem A.1: Plot of $e^{-t}$','Interpreter','latex');
xlabel('t');
ylabel('$e^{-t}$','Interpreter','latex');
grid on;

f_exp_pos = @(t) exp(t);
subplot(2, 1, 2);
plot(t1, f_exp_pos(t1));
title('Problem A.1: Plot of $e^{t}$','Interpreter','latex');
xlabel('t');
ylabel('$e^{t}$','Interpreter','latex');
grid on;

% Problem A.2 & A.3
t2 = -2:1:2;
y2 = f_exp(t2);
figure('Name','Problem A.2 & A.3');
stem(t2, y2, 'filled');
hold on;
plot(t1, f_exp(t1), 'r--');
title('Problem A.2 & A.3: Comparison of Plots');
xlabel('t');
ylabel('$e^{-t}$','Interpreter','latex');
legend('Sampled Points (t = -2:2)', 'Continuous Function');
grid on;
disp('Problem A.3: The plot from A.2 shows 5 discrete points, while the plot from A.1 is a smooth curve due to a much denser time vector.');
disp(' ');

%% =====================================================================
% B. Time Shifting and Time Scaling
% ======================================================================
disp('Executing Section B: Time Shifting & Scaling...');

% Problem B.1
p = @(t) (t >= 0 & t < 2);
t_b = -1:0.01:5;
figure('Name', 'Section B Plots');
subplot(2, 2, 1);
plot(t_b, p(t_b), 'LineWidth', 2);
title('Problem B.1: p(t)');
xlabel('t');
ylim([-0.2, 1.2]);
grid on;

% Problem B.2
r = @(t) t .* p(t);
n = @(t) r(t) + r(-t + 2);
subplot(2, 2, 2);
plot(t_b, r(t_b), 'LineWidth', 2);
title('Problem B.2: r(t) = t*p(t)');
xlabel('t');
grid on;
subplot(2, 2, 3);
plot(t_b, n(t_b), 'LineWidth', 2);
title('Problem B.2: n(t) = r(t) + r(-t+2)');
xlabel('t');
grid on;

% Problem B.3, B.4, & B.5
n1 = @(t) n(0.5 * t);
n2 = @(t) n1(t + 0.5); 
n3 = @(t) n(t + 0.25);
n4 = @(t) n3(0.5 * t);
figure('Name', 'Problem B.3-B.5 Plots');
subplot(2, 1, 1);
plot(t_b, n2(t_b));
title('Problem B.3: $n_2(t) = n(0.5(t+0.5))$','Interpreter','latex');
xlabel('t');
grid on;
subplot(2, 1, 2);
plot(t_b, n4(t_b));
title('Problem B.4: $n_4(t) = n(0.5t+0.25)$','Interpreter','latex');
xlabel('t');
grid on;
disp('Problem B.5: The plots for n2(t) and n4(t) are identical because both simplify to n(0.5*t + 0.25).');
disp(' ');

%% =====================================================================
% C. Visualizing Operations and Algorithm Vectorization
% ======================================================================
disp('Executing Section C: Visualizing Operations & Vectorization...');

% Problem C.1 & C.2
f = @(t) exp(-2*t) .* cos(4*pi*t);
u = @(t) (t >= 0);
g = @(t) f(t) .* u(t);
t_c1 = -2:0.01:4;
figure('Name', 'Problem C.1 & C.2');
subplot(2,1,1);
plot(t_c1, g(t_c1));
title('Problem C.1: $g(t) = e^{-2t}\cos(4\pi t)u(t)$','Interpreter','latex');
xlabel('t');
grid on;

s = @(t) g(t + 1);
subplot(2,1,2);
plot(t_c1, s(t_c1));
title('Problem C.2: $s(t) = g(t+1)$','Interpreter','latex');
xlabel('t');
grid on;

% Problem C.3 & C.4
t_c3 = 0:0.01:4;
alpha = [1, 3, 5, 7];
t_matrix = repmat(t_c3, length(alpha), 1);
alpha_matrix = repmat(alpha', 1, length(t_c3));
s_alpha_matrix = exp(-2) * exp(-alpha_matrix .* t_matrix) .* cos(4*pi*t_matrix);
figure('Name', 'Problem C.3');
plot(t_c3, s_alpha_matrix);
title('Problem C.3: $s_{\alpha}(t) = e^{-2}e^{-\alpha t}\cos(4\pi t)u(t)$','Interpreter','latex');
xlabel('t');
ylabel('$s_{\alpha}(t)$','Interpreter','latex');
legend('$\alpha = 1$', '$\alpha = 3$', '$\alpha = 5$', '$\alpha = 7$','Interpreter','latex');
grid on;
matrix_size = size(s_alpha_matrix);
disp('Problem C.4:');
fprintf('The size of the matrix generated in C.3 is %d x %d.\n', matrix_size(1), matrix_size(2));
disp(' ');

%% =====================================================================
% D. Array Indexing
% ======================================================================
disp('Executing Section D: Array Indexing...');

% Load data from .mat file
load('ELE532_Lab1_Data.mat');
disp("'ELE532_Lab1_Data.mat' file loaded successfully.");
disp(' ');

% --- Problem D.1 ---
disp('--- Problem D.1 ---');
A = [ 0.5377 -1.3077 -1.3499 -0.2050;
      1.8339 -0.4336  3.0349 -0.1241;
     -2.2588  0.8622  3.5784 -0.0631;
      0.3426  0.7254 -0.0631  1.4172;
      1.4090 -0.4447  0.7147  1.4897 ];
disp('Original Matrix A:');
disp(A);

disp('(a) A(:) reshapes A into a single column vector.');
disp(A(:));

disp('(b) A([2 4 7]) extracts the 2nd, 4th, and 7th elements using linear indexing.');
disp(A([2 4 7]));

disp('(c) A >= 0.2 returns a logical matrix of the same size as A.');
disp(A >= 0.2);

disp('(d) A(A >= 0.2) extracts elements from A where the condition is true.');
disp(A(A >= 0.2));

A_temp = A; % Create a temporary copy for this operation
A_temp(A_temp >= 0.2) = 0;
disp('(e) A(A >= 0.2) = 0 sets elements meeting the condition to zero.');
disp(A_temp);
disp(' ');

% --- Problem D.2 ---
disp('--- Problem D.2 ---');
B_loop = B; 
tic;
for j = 1:size(B_loop, 2)
    for i = 1:size(B_loop, 1)
        if abs(B_loop(i, j)) < 0.01
            B_loop(i, j) = 0;
        end
    end
end
time_loop = toc;
fprintf('Nested for-loops execution time: %f seconds\n', time_loop);

B_indexed = B;
tic;
B_indexed(abs(B_indexed) < 0.01) = 0;
time_indexed = toc;
fprintf('Indexing execution time: %f seconds\n', time_indexed);
fprintf('Comparison: Indexing is approximately %.2f times faster.\n', time_loop / time_indexed);
disp(' ');

% --- Problem D.3 ---
disp('--- Problem D.3 ---');
disp('Executing simple audio compression algorithm...');
fs = 8000;
xaudio_processed = x_audio;
threshold = 0.05;

below_threshold_indices = abs(xaudio_processed) < threshold;
num_zeros = sum(below_threshold_indices);
compression_ratio = (num_zeros / length(x_audio)) * 100;
xaudio_processed(below_threshold_indices) = 0;

fprintf('Threshold used: %.2f\n', threshold);
fprintf('Number of samples set to zero: %d / %d\n', num_zeros, length(x_audio));
fprintf('Compression ratio (percentage of zeroed samples): %.2f%%\n', compression_ratio);

% To listen to the audio, uncomment the lines below
% disp('Playing original audio...');
% sound(x_audio, fs);
% pause(3);
% disp('Playing processed audio...');
% sound(xaudio_processed, fs);

disp(' ');
disp('All tasks completed.');