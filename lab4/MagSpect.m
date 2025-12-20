function [ output_args ] = MagSpect( x )
% MAGSPECT ... Utility function to simplify plotting the magnitude spectrum.
%
% This is the MODIFIED version that does NOT require the Signal Processing Toolbox.

%==========================================================================
% Default values
%==========================================================================
Fs   = 32000;   % default sampling frequency
Nfft = 1024;    % default FFT size

%==========================================================================
% Set up the frequency vector
%==========================================================================
ff = [ -(Nfft/2) : 1 : (Nfft/2)-1 ] * (Fs/Nfft);

%==========================================================================
% Compute the spectrum of x(t) using Nfft-point FFT 
%==========================================================================
Xspect = fft( x, Nfft );

%==========================================================================
% plot the magnitude spectrum 
%==========================================================================

% --- MODIFICATION ---
% Original line (line 36): plot( ff, db( abs( fftshift(Xspect) ) ) );
% The 'db' function requires the Signal Processing Toolbox.
% We replace it with the standard formula: 20*log10(Magnitude)

% Calculate magnitude and add a small value (eps) to avoid log(0) = -Inf
magnitude_spectrum = abs( fftshift(Xspect) ) + eps; 

plot( ff, 20*log10( magnitude_spectrum ) );
% --- END MODIFICATION ---

set(gca,'XLim',[-Fs/2 Fs/2]);
    xlabel('Frequency [Hz]'); ...
    ylabel('Magnitude [dB]'); ...
    grid on; 

end