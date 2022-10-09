function get_plot(data, x, y, cmap, xlab, ylab, ti)
%GET_PLOT Plot imagesc with specified parameters
%   Plots data with specified colormap, x-label, y-label and title.

imagesc(x, y, data);    % plots data
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
