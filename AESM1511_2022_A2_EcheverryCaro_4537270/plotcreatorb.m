function plotcreatorb(data,cmap, xlab, ylab, ti, xvalue,yvalue)
% plotcreator plots imagesc with specified colormap, x-label, y-label and title.
 
imagesc(data,'xData',xvalue,'yData',yvalue);    % plots data
set(gca,'linewidth',2); % sets linewidth
colormap(cmap);         % set colormap to parameter
colorbar;               % plot colorbar
t = title(ti);          % set title to parameter
lx = xlabel(xlab);      % set xlabel to parameter
ly = ylabel(ylab);      % set ylabel to parameter
t.FontSize = 14;        % set font sizes
lx.FontSize = 14;
ly.FontSize = 14;
end