% ELE 532 Signals & Systems I - Laboratory Assignment 3

% --- Environment Initialization ---
clear all;
close all;
clc;

%% ========================================================================
%  Create Main Figure and Tab Group
% =========================================================================

mainFig = figure('Name', 'ELE 532 Lab 3 Results - Fourier Series Analysis', 'NumberTitle', 'off', 'WindowState', 'maximized');
tabGroup = uitabgroup('Parent', mainFig);
fprintf('Main figure window created.\n\n');

%% ========================================================================
%  Problem A.4: Generate and Plot Spectra (Magnitude and Phase)
% =========================================================================

fprintf('Executing Problem A.4: Generating and Plotting Spectra...\n');

n_ranges = {(-5:5), (-20:20), (-50:50), (-500:500)};
signal_types = {'x1', 'x2', 'x3'};

for i = 1:length(signal_types)
    signal_name = signal_types{i};
    tab_spec = uitab('Parent', tabGroup, 'Title', [signal_name ' Spectra']);
    
    for j = 1:length(n_ranges)
        n_vec = n_ranges{j};
        Dn = calculate_Dn(n_vec, signal_name);
        
        % Subplot for Magnitude (Left Column)
        subplot(4, 2, 2*j-1, 'Parent', tab_spec); % Position: 1, 3, 5, 7
        stem(n_vec, abs(Dn), 'filled');
        title(['Magnitude |D_n| for N = ' num2str(max(n_vec))]);
        xlabel('n');
        ylabel('Magnitude');
        grid on;

        % Subplot for Phase (Right Column)
        subplot(4, 2, 2*j, 'Parent', tab_spec); % Position: 2, 4, 6, 8
        stem(n_vec, angle(Dn), 'filled', 'r'); % angle() calculates phase in radians
        title(['Phase ∠D_n for N = ' num2str(max(n_vec))]);
        xlabel('n');
        ylabel('Phase (radians)');
        ylim([-pi*1.1, pi*1.1]); % Set y-axis limits for phase
        grid on;
    end
end
fprintf('All magnitude and phase spectra have been plotted.\n\n');

%% ========================================================================
%  Problem A.6: Reconstruct Signals from Fourier Series (This part remains the same)
% =========================================================================
fprintf('Executing Problem A.6: Reconstructing Signals...\n');
t = -300:1:300;
for i = 1:length(signal_types)
    signal_name = signal_types{i};
    tab_recon = uitab('Parent', tabGroup, 'Title', [signal_name ' Reconstruction']);
    ax_recon = axes('Parent', tab_recon);
    hold(ax_recon, 'on');
    legend_entries = {};
    for j = 1:length(n_ranges)
        n_vec = n_ranges{j};
        Dn = calculate_Dn(n_vec, signal_name);
        x_reconstructed = reconstruct_signal(Dn, n_vec, signal_name, t);
        plot(ax_recon, t, real(x_reconstructed), 'LineWidth', 1.5);
        legend_entries{j} = ['N = ' num2str(max(n_vec))];
    end
    hold(ax_recon, 'off');
    title(ax_recon, ['Reconstructed Signal ' signal_name '(t) with Varying Number of Coefficients']);
    xlabel(ax_recon, 'Time (t)');
    ylabel(ax_recon, 'Amplitude');
    legend(ax_recon, legend_entries);
    grid(ax_recon, 'on');
    xlim(ax_recon, [-50, 50]);
end
fprintf('All signals have been reconstructed and plotted.\n');
fprintf('\nLab 3 script execution complete.\n');

%% ========================================================================
%  Helper Functions (These remain the same)
% =========================================================================
function Dn = calculate_Dn(n_range, signal_type)
    Dn = zeros(size(n_range));
    switch signal_type
        case 'x1'
            Dn(n_range == 1) = 0.25; Dn(n_range == -1) = 0.25;
            Dn(n_range == 3) = 0.5; Dn(n_range == -3) = 0.5;
        case 'x2'
            Dn(n_range == 0) = 0.5;
            n_nonzero = n_range(n_range ~= 0);
            Dn(n_range ~= 0) = 0.5 * sin(pi * n_nonzero / 2) ./ (pi * n_nonzero / 2);
        case 'x3'
            Dn(n_range == 0) = 0.4;
            n_nonzero = n_range(n_range ~= 0);
            Dn(n_range ~= 0) = 0.4 * sin(pi * 2 * n_nonzero / 5) ./ (pi * 2 * n_nonzero / 5);
    end
end
function x_rec = reconstruct_signal(Dn, n_range, signal_type, t)
    if strcmp(signal_type, 'x1') || strcmp(signal_type, 'x2')
        w0 = pi/10;
    else, w0 = 2*pi/25;
    end
    x_rec = zeros(size(t));
    for i = 1:length(n_range)
        n = n_range(i); D_n = Dn(i);
        x_rec = x_rec + D_n * exp(1j * n * w0 * t);
    end
end