% function [frq, amp, phase] = simpleFFT( signal, ScanRate)
% Purpose: perform an FFT of a real-valued input signal, and generate the single-sided 

% output, in amplitude and phase, scaled to the same units as the input.

%inputs: 

%    signal: the signal to transform

%    ScanRate: the sampling frequency (in Hertz)

% outputs:

%    frq: a vector of frequency points (in Hertz)

%    amp: a vector of amplitudes (same units as the input signal)

%    phase: a vector of phases (in radians)

function [frq, real_part, imag_part, amp, phase] = simpleFFT( signal, ScanRate)

n = length(signal); 

z = fft(signal, n); %do the actual work

%generate the vector of frequencies

halfn = floor(n / 2)+1;

deltaf = 1 / ( n / ScanRate);

%frq = (-halfn+1:(halfn-1)) * deltaf;
frq = (0:(halfn-1)) * deltaf;

% convert from 2 sided spectrum to 1 sided

%(assuming that the input is a real signal)


amp(1) = abs(z(1)) ./ (n);

amp(2:(halfn-1)) = abs(z(2:(halfn-1))) ./ (n / 2); 

amp(halfn) = abs(z(halfn)) ./ (n); 

% amp = [fliplr(amp(2:end)), amp];

phase = angle(z(1:halfn));

% phase = [fliplr(phase(2:end)), phase];

[real_part, imag_part] = pol2cart(phase, amp);

end