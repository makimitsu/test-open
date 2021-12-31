function [] = plot_B_z_in_time(B_z,ch_dist,t_start,t_end)
% plot B_z in time space
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   2d array of int: ch_dist, keep track of channel locations
%   int: t_start
%   int: t_end
r = 8;
col1 = 14;
col2 = 14;

figure('Position', [10 10 1200 900])
for i=1:r
    for j=1:col1
        subplot(r,col1,(i-1)*col1+j)
        plot(reshape(B_z(i,j,t_start:t_end),[],1))
        title(num2str(ch_dist(i,j)));
        ylim([-0.04 0.04]);
    end
end

figure('Position', [10 10 1200 900])
for i=1:r
    for j=col1+1:col1+col2
        subplot(r,col2,(i-1)*col2+j-col1)
        plot(reshape(B_z(i,j,t_start:t_end),[],1))
        title(num2str(ch_dist(i,j)));
        ylim([-0.04 0.04]);
    end
end
end