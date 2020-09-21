function [x_interp,y_interp,interp_data] = PiecewiseLinearInterpolate(y_function,x_interval,N)

    dx = 1/N*(x_interval(2) - x_interval(1));
    
    x_interp = x_interval(1):dx:x_interval(2); % length N + 1
    y_interp = y_function(x_interp); % length N + 1
    
    interp_data = cell(1,N);
    
    for i = 1:N
        slope = (y_interp(i+1) - y_interp(i))/(x_interp(i+1) - x_interp(i));
        intercept = y_interp(i) - slope*x_interp(i);
        interp_data{i}.slope = slope;
        interp_data{i}.intercept = intercept;
        interp_data{i}.xstart = x_interp(i);
        interp_data{i}.ystart = y_interp(i);
    end

end

