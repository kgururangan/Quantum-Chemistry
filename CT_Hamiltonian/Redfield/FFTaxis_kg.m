function [y] = FFTaxis_kg(t,centered,linear,real)
    % Constructs the frequency axis corresponding to time axis t
    % Examples:
    % 1-sided half frequency axis for real signals in Hz - f = FFTaxis_kg(t,0,1,1)
    % 1-sided full frequency axis for complex signals in Hz - f = FFTaxis_kg(t,0,1,0)
    % 2-sided frequency axis for complex signals in rad/s - omega =
    % FFTaxis_kg(t,1,0,0) (MUST fftshift() spectrum to match this axis!!!)


    dt = abs(t(2)-t(1));
    fs = 1/dt;

    if centered == 1
        % MUST fftshift() spectrum to match with centered axis!!!
        disp('Remember to fftshift() spectrum to match with centered axis!')
        N = [floor(-length(t)/2):floor(length(t)/2)-1];
    else
        N = [0:1:length(t)-1];
    end
    
    if centered ~= 1 && real == 1
        N = [0:1:length(t)/2-1];
    end

    if linear == 1 % Hz
        % e.g. cos(2*pi*f0*t) -> peak appears at f0
        y = fs/(length(t)-1)*N;
    else % rad/s
        % e.g. cos(f0*t) -> peak appears at f0
        y = 2*pi*fs/(length(t)-1)*N;
    end

end