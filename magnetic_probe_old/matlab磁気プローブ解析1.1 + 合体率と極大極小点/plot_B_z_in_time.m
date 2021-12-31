function [] = plot_B_z_in_time(B_z,ch_dist,t_start,t_end)
% plot B_z in time space
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   2d array of int: ch_dist, keep track of channel locations
%   int: t_start
%   int: t_end

figure('Position', [0 0 1500 1500])
sgtitle('B(T) vs t(us)')
for i = 1:7
    for j = 1:13
        subplot(7,13,(i-1)*13+j)
        plot(squeeze(B_z(i,j,:)));
        xlim([t_start t_end])
        ylim([-0.02 0.1])
        title(num2str(ch_dist(i,j)));
    end
end
sgtitle('B(T) vs t(us)')
figure('Position', [0 0 1500 1500])
for i = 1:7
    for j = 1:12
        subplot(7,12,(i-1)*12+j)
        plot(squeeze(B_z(i,j+13,:)));
        xlim([t_start t_end])
        ylim([-0.02 0.1])
        xlabel('t (us)');
        ylabel('B (T)');
        title(num2str(ch_dist(i,j+13)));
    end
end
end