function scrollImagesc_plots(x, y, data, clim)

%cmap = jet_white;
cmap = parula;

% clim = 'manual' or 'auto'
% 'manual' = colorbar set at range of first image
% 'auto' = colorbar adjusts for current image

numImages = size(data,3);
	
	switch nargin
		case 1
			data = x;
			x = 1:size(data, 2);
			y = 1:size(data, 1);
			clim = 'manual';
		case 3
			clim = 'manual';
		case 4
			
		otherwise
			InvalidNumInputArgs('scrollImagesc');
	end
	figure('WindowScrollWheelFcn', @wheelCallback, ...
		'WindowKeyPressFcn', @keyCallback);
	a = 1;
        contour(x,y,data(:,:,a),20,'Linewidth',2); hold off;
        xlabel('-\omega_3/cm^{-1}','FontSize',16)
        ylabel('\omega_1/cm^{-1}','FontSize',16)
        set(gca,'Ydir','normal','FontSize',16,'Linewidth',1.5,'Box','off','xtick',[1.2e4,1.25e4,1.3e4],'ytick',[1.2e4,1.25e4,1.3e4])
        colorbar
        colormap(cmap)
        %colormap(brewermap([],'RdBu'))
        %colormap(redblue_alt)
        axis square
        axis([1.2e4,1.3e4, 1.2e4,1.3e4])
        grid on
    title(sprintf('%u / %u',a,numImages));
	
	function wheelCallback(src, cbd)
		a = a+cbd.VerticalScrollCount;
		if a<1
			a=1;
		elseif a>size(data, 3)
			a=size(data, 3);
		end
		hold all;
		caxis(clim);
		cla;
		%imagesc(x, y, data(:,:,a));
        contour(x,y,data(:,:,a),20,'Linewidth',2); hold off;
        xlabel('-\omega_3/cm^{-1}','FontSize',16)
        ylabel('\omega_1/cm^{-1}','FontSize',16)
        set(gca,'Ydir','normal','FontSize',16,'Linewidth',1.5,'Box','off','xtick',[1.2e4,1.25e4,1.3e4],'ytick',[1.2e4,1.25e4,1.3e4])
        colorbar
        colormap(cmap)
        axis square
        axis([1.2e4,1.3e4, 1.2e4,1.3e4])
        grid on
        title(sprintf('%u/%u',a,numImages));
	end

	function keyCallback(src, cbd)
		if double(cbd.Character) == 30
			a = a+1;
		elseif double(cbd.Character) == 31
			a = a-1;
		end
		if a<1
			a=1;
		elseif a>size(data, 3)
			a=size(data, 3);
		end
		hold all;
		caxis(clim);
		cla;
        contour(x,y,data(:,:,a),20,'Linewidth',2); hold off;
        xlabel('-\omega_3/cm^{-1}','FontSize',16)
        ylabel('\omega_1/cm^{-1}','FontSize',16)
        set(gca,'Ydir','normal','FontSize',16,'Linewidth',1.5,'Box','off','xtick',[1.2e4,1.25e4,1.3e4],'ytick',[1.2e4,1.25e4,1.3e4])
        colorbar
        colormap(cmap)
        axis square
        axis([1.2e4,1.3e4, 1.2e4,1.3e4])
        grid on
        title(sprintf('%u/%u',a,numImages));
	end

end

