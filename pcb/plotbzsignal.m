function plotbzsignal(y_upper_lim, col2, col1, t_end, p_ch, y_lower_lim, t_start, bz, ok, r_ch, r)
figure('Position', [10 10 1200 900])
%bz_s=smoothdata(bz,1);
for i=1:r
    for j=1:col1
        subplot(r,col1,(i-1)*col1+j)
        if ok(r_ch*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j))
            %plot(t_start:t_end,bz_s(t_start:t_end,r_ch*(i-1)+j),'r')
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j),'r:')
        end   
        title(num2str(p_ch(i,j)));
        %xlim([t_start t_end]);
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
        %ylim([-0.02 0.04]);
    end
end

figure('Position', [10 10 1200 900])
for i=1:r
    for j=col1+1:col1+col2
        subplot(r,col2,(i-1)*col2+j-col1)
        if ok(r_ch*(i-1)+j)==1 
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j))
            %plot(t_start:t_end,bz_s(t_start:t_end,r_ch*(i-1)+j),'r')
        else 
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j),'r:')
        end   
        title(num2str(p_ch(i,j)));
        %xlim([t_start t_end]);
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
        %ylim([-0.02 0.04]);
    end
end
end

%パラメータの説明
% r = 5;%プローブ本数＝グラフ出力時の縦に並べる個数
% col1 = 12;%1枚目のグラフ出力時の横に並べる個数
% col2 = 13;%2枚目のグラフ出力時の横に並べる個数
% y_upper_lim = 0.1;%縦軸プロット領域（b_z上限）
% y_lower_lim = -0.1;（b_z下限）
% t_start=455;%横軸プロット領域（開始時間）
% t_end=520;%横軸プロット領域（終了時間）
% r_ch=col1+col2;%r方向から挿入した各プローブのチャンネル数


