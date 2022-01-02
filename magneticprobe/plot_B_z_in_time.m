function [] = plot_B_z_in_time(B_z,ch_dist,t_start,t_end)
% plot B_z in time space
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   2d array of int: ch_dist, keep track of channel locations
%   int: t_start
%   int: t_end
r = size(B_z,1);
col1 = round(size(B_z,2)/2);
col2 = size(B_z,2)-col1;
y_upper_lim = 0.08;
y_lower_lim = -0.08;

figure('Position', [10 10 1200 900])
for i=1:r
    for j=1:col1
        subplot(r,col1,(i-1)*col1+j)
        plot(reshape(B_z(i,j,t_start:t_end),[],1))
        title(num2str(ch_dist(i,j)));
        ylim([y_lower_lim y_upper_lim]);
    end
end

figure('Position', [10 10 1200 900])
for i=1:r
    for j=col1+1:col1+col2
        subplot(r,col2,(i-1)*col2+j-col1)
        plot(reshape(B_z(i,j,t_start:t_end),[],1))
        title(num2str(ch_dist(i,j)));
        ylim([y_lower_lim y_upper_lim]);
    end
end
end