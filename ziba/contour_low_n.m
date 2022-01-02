function contour_low_n(low_n_signal,t_start,t_end)

% ## set y axis ##
not_rad_x = [8.6, 68, 113, 143, 188, 233, 263, 323];
x = not_rad_x*2*pi/360;

% ## set x axis ##
t = linspace(t_start,t_end,t_end*10-t_start*10+1);

% ## plot countour figure ##
figure('Position', [10 10 800 400]);
contourf(t,x,low_n_signal(t_start*10:t_end*10,:)','LineStyle','none');
xlabel('time[us]','Fontsize',12,'FontWeight','bold');
ylabel('toroidal position[rad]','Fontsize',12,'FontWeight','bold');
yticks(x);
c = colorbar;
colormap turbo;
%colorbar('Fontsize',11,'FontWeight','bold')
c.Label.String = 'B[T]';
c.Label.FontSize = 12;
c.Label.FontWeight = 'bold';
%caxis([-0.15 0.15]);
