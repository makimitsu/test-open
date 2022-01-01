function [] = plot_fitrate(B_z,r_probe,shot)
% plot fit rate(not sure correct words) and magnetic surface
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z

time = 460:500;
x = 1:550;
y = 1:550;

hold on
                                
for i = time
    psi = get_psi(B_z,r_probe,i);
    max_psi = max(psi,[],1);
    plot(max_psi);
    TF = islocalmin(max_psi);
    if isempty(max_psi(TF)) == 1
        psi_com = 1;
    else 
        psi_com = min(max_psi(TF));  %psi_com: psi_common(after combine)
    end
    for k = 1:length(max_psi)
        if psi_com == max_psi(k) 
            if TF(k) == 1
                j = k;
            end
        end
    end
    plot(j, max_psi(j),'r*')
    if j == 0
        y(i) = 0;
    else
        [psi_pri, k] = max(max_psi(1:j));
        [psi_pri_lat, l] = max(max_psi(j:25));
        if psi_pri > psi_pri_lat
            psi_pri = psi_pri_lat;
            k = l;
        end
        if psi_com/psi_pri > 1
            y(i) = 1;
        elseif psi_pri == 0
            y(i) = 1;
        else
            y(i) = psi_com/psi_pri;
        end
    end
    plot(k, max_psi(k), '*b');
end

hold off

xlabel('r');
ylabel('max_psi(Wb)');

f = figure('name',['shot', num2str(shot)]);

plot(x, y, '-o');
title('');
xlabel('time(us)');
ylabel('fitrate');
axis([460 500 0 1.0]);
end