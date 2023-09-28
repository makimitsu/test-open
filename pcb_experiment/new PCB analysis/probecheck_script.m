%%%%%%%%%%%%%%%%%%%%%%%%
%200ch用新規pcbプローブと装置の磁場信号極性チェック
%dtacqのshot番号を直接指定する場合
%%%%%%%%%%%%%%%%%%%%%%%%

figure_switch = ["on","on","off","off","off","off"];%bz1,bz2,bt1,bt2,bz_vs_z,bt_vs_z

r = 7;%プローブ本数＝グラフ出力時の縦に並べる個数
col = 10;%グラフ出力時の横に並べる個数
y_upper_lim = 0.03;%縦軸プロット領域（b_z上限）
y_lower_lim = -0.03;%縦軸プロット領域（b_z下限）
t_start=300;%横軸プロット領域（開始時間）
t_end=800;%横軸プロット領域（終了時間）
t_start=1;
t_end=1000;
t = 494;
% z_probe_pcb = [-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17];
z_probe_pcb = [-0.2975 -0.255 -0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17 0.255 0.2975];
n_z = length(z_probe_pcb);
% r_ch=col1+col2;%r方向から挿入した各プローブのチャンネル数

f1=figure(Visible=figure_switch(1));
f1.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bz(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i-1)+j)-1));
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bz signal probe1-5')


f2=figure(Visible=figure_switch(2));
f2.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bz(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i+r-1)+j)-1));
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bz signal probe6-10')

f3=figure(Visible=figure_switch(3));
f3.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bt(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i-1)+j)));
        xticks([t_start t_end]);
        %ylim([-0.2 0.2]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bt signal probe1-5')

f4=figure(Visible=figure_switch(4));
f4.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bt(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i+r-1)+j)));
        xticks([t_start t_end]);
        %ylim([-0.2 0.2]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bt signal probe6-10')

% saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_02','.png'))
% close

% 横軸z, 縦軸Bzのプロット
f5=figure(Visible=figure_switch(5));
f5.WindowState = 'maximized';
styles = ["-*r","-^r","-vr","-<b","->b","-ob","-squarek","-diamondk","-pentagramk","-hexagramk"];
tiles = tiledlayout(2,1);
sgtitle(strcat('t=',num2str(t),' us'))
nexttile
hold on
for i=1:10
    %14本のプローブの一番内側の点：zline == [1,11,21,31,41,51,61,71,81,91,101,111,121,131]
    zline=(1:10:n_z*10-9)+(i-1);%1~10,11~20....
    % r = 0.21mの信号,中心6本
    %zline = [51,71,91,111,131,151];
    bz_zline=bz(t,zline);
    bz_zline(ok_bz(zline)==false)=NaN;
    plot(z_probe_pcb,bz_zline,styles(i),'Color','k','MarkerSize',12)
    clear bz_zline
end
hold off
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
title('Before interpolation')

Bz_interped = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
nexttile
hold on
for i=1:10
    zline=(1:10:n_z*10-9)+(i-1);
    bz_zline=Bz_interped(zline);
    plot(linspace(min(zpos_bz),max(zpos_bz),n_z),bz_zline,styles(i),'Color','k','MarkerSize',12)
    clear bz_zline
end
hold off
xlabel('z [m]')
ylabel('Bz')
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
title('After interpolation')
tiles.TileSpacing = 'compact';
tiles.Padding = 'compact';

% 横軸z, 縦軸Btのプロット
f6=figure(Visible=figure_switch(6));
f6.WindowState = 'maximized';
styles = ["-*","-^","-v","-<","->","-o","-square","-diamond","-pentagram","-hexagram"];
tiles = tiledlayout(2,1);
sgtitle(strcat('t=',num2str(t),' us'))
nexttile
hold on
for i=1:10
    zline=(1:10:n_z*10-9)+(i-1);
    bt_zline=bt(t,zline);
    bt_zline(ok_bt(zline)==false)=NaN;
    plot(z_probe_pcb,bt_zline,styles(i),'Color','k','MarkerSize',12)
    clear bz_zline
end
hold off
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
title('Before interpolation')

Bt_interped = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, bt, ok_bt, t);
nexttile
hold on
for i=1:10
    zline=(1:10:n_z*10-9)+(i-1);
    bt_zline=Bt_interped(zline);
    plot(linspace(min(zpos_bz),max(zpos_bz),n_z),bt_zline,styles(i),'Color','k','MarkerSize',12)
    clear bz_zline
end
hold off
xlabel('z [m]')
ylabel('Bt')
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
title('After interpolation')
tiles.TileSpacing = 'compact';
tiles.Padding = 'compact';

hidden = find(figure_switch == 'off');
figures = [f1,f2,f3,f4,f5,f6];
for i = hidden
    close(figures(i));
end

