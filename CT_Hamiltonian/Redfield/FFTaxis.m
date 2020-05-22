function [ y, x, N ] = FFTaxis( n, dx, centered )
%FFTaxis Calculates the conjugate axis for the FFT of a data series
%	For function DATA(x), its fast Fourier transform is denoted DATA'(y)
%	where DATA, DATA', x, and y are all of length n.
%
%	If x represents time, then y represents frequency or vice versa.
%
%   When going from time to frequency, the resulting axis is in units of
%   linear frequency.
%
%	y = FFTaxis(x) returns the axis y conjugate to axis x along with
%		unit-step axis N given a vector x.  If x is zero-centered, so will
%		be y.
%
%	y = FFTaxis(n, dx) returns the axis y conjugate to axis x given the
%		scalar axis length n and a scalar x-axis step dx.
%
%	y = FFTaxis(x, CENTERED) returns the axis y conjugate to axis x given
%		a vector x.  Boolean value CENTERED determines if the returned y
%		axis is zero-centered or not.
%
%	y = FFTaxis(DATA, x) returns the axis y conjugate to axis x given
%		a DATA vector and an axis vector x, both of length n.  If x is
%		zero-centered, so will be y.
%
%	y = FFTaxis(DATA, dx) returns the axis y conjugate to axis x given
%		a DATA vector (of length n) and an scalar axis step dx.
%
%	y = FFTaxis(n, x) returns the axis y conjugate to axis x given
%		a scalar axis length n and a vector axis x.  If x is zero-centered,
%		so will be y.
%
%	y = FFTaxis(n, dx, CENTERED) returns the axis y conjugate to axis x
%		given the scalar axis length n and a scalar x-axis step dx.
%		Boolean value CENTERED determines if the returned y axis is
%		zero-centered or not.
%
%	[y, x, N] = FFTaxis(...) returns the axis y conjugate to axis x along
%		with unit-step axis N.
%
%	See also FFTpad, fft, ifft, fft2, ifft2, fftn, ifftn.
	
	switch nargin
		case 1
			% FFTaxis(x)
			dx = n;
			centered = dx(floor(length(dx)/2+1))==0;
		case 2
			if ~isscalar(n) && islogical(dx) && isscalar(dx)
				% FFTaxis(x, centered)
				centered = dx;
				dx = n;
			elseif ~isscalar(n) && ~isscalar(dx)
				% FFTaxis(data, x)
				centered = dx(floor(length(dx)/2+1))==0;
			elseif ~isscalar(n)
				% FFTaxis(data, dx)
				centered = true;
			elseif isscalar(n) && ~isscalar(dx)
				% FFTaxis(n, x)
				centered = dx(floor(length(dx)/2+1))==0;
			else
				centered = true;
			end
	end
	
	if ~isscalar(n)
		n = length(n);
	end
	
	if ~isscalar(dx)
		dx = mean(diff(dx(:)));
	end
	
	if centered
		N = floor(-n/2):1:floor(n/2 - 1);
	else
		N = 0:(n-1);
	end
	
	if dx<0
		N = -flip(N);
	end
	
	y = N/(dx*n);
	x = N*dx;

end

% function [ y, x, N ] = FFTaxis( n, dx, centered )
% %FFTaxis Calculates the conjugate axis for the FFT of a data series
% %   When going from time to frequency, the resulting axis is in units of
% %   linear frequency.
% 
% 	switch nargin
% 		case 1
% 			% n is the pre-FFT axis
% 			dx = n;
% 			centered = dx(length(dx)/2+1)==0;
% 		case 2
% 			centered = true;
% 	end
% 	
% 	if ~isscalar(n)
% 		n = length(n);
% 	end
% 	
% 	if ~isscalar(dx)
% 		dx = mean(diff(dx(:)));
% 	end
% 	
% 	if centered
% 		N = floor(-n/2):1:floor(n/2 - 1);
% 	else
% 		N = 0:(n-1);
% 	end
% 	
% 	y = N/(dx*n);
% 	x = N*dx;
% 
% end

% function [ f, t, N ] = FFTaxis( n, dt, centered )
% %FFTaxis Calculates the conjugate axis for the FFT of a data series
% 
% 	switch nargin
% 		case 2
% 			centered = true;
% 	end
% 	
% 	if ~isscalar(n)
% 		n = length(n);
% 	end
% 	
% 	if centered
% 		N = floor(-n/2):1:floor(n/2 - 1);
% 	else
% 		N = 0:(n-1);
% 	end
% 	
% 	f = N/(dt*n);
% 	t = N*dt;
% 
% end

%==================ORIGINAL CODE======================%
% function [ freqAxis ] = FFTaxis( numPoints, timeStep )
% %FFTaxis Calculates the angular frequency axis for the FFT of a time series
% %   Detailed explanation goes here
% 
%     samplingFreq = (2*pi)/timeStep;
%     freqStep = samplingFreq/numPoints;
% %     indices = -(numPoints/2 - 1):numPoints/2;
% 	indices = -(numPoints/2):(numPoints/2 - 1);
%     freqAxis = indices*freqStep;
% 
% %     freqAxis = zeros(numPoints,1);
% %     fs = 1/timeStep;
% %     fn = fs/2;
% %     freqStep = fs/(numPoints-1);
% %     for i = 1:numPoints
% %         freqAxis(i) = freqStep*(i-1);
% %     end
% end
