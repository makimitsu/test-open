%%%%%%%%%%%%%%%%%%%%%%%%
%125ch用 akimitsu pcbプローブと装置の磁場信号極性チェック
%dtacqのshot番号を直接指定する場合
%%%%%%%%%%%%%%%%%%%%%%%%
pathname.fig = "C:\Users\uswk0\OneDrive - The University of Tokyo\data\figure\"; % plot画像の保存先

figure_switch = ["on","on","on","on"];%bz1,bz2,bz_vs_z, bz_vs_r

r = 5;%プローブ本数＝グラフ出力時の縦に並べる個数
col1 = 13;%グラフ出力時の横に並べる個数 1枚目
col2 = 12;%グラフ出力時の横に並べる個数 2枚目
r_ch = col1 + col2;% 2枚のグラフの合計（１本のプローブのch数）
y_upper_lim = 0.05;%縦軸プロット領域（b_z上限）
y_lower_lim = -0.05;%0.03;%縦軸プロット領域（b_z下限）
t_start=300;%横軸プロット領域（開始時間）
t_end=800;%横軸プロット領域（終了時間）
r_shift = 0.1;
t_start=1;
t_end=1000;

t = 478;
z_probe_pcb    = [-0.0525 -0.021 0 0.021 0.0525]; % ch1-25, ch26-50, ch51-76, ch77-101, ch102-126

r_probe_pcb    = [0.0600 0.0800 0.1000 0.1200 0.1400 0.1500...
                  0.1550 0.1600 0.1650 0.1700 0.1750 0.1800...
                  0.1850 0.1900 0.1950 0.2000 0.2050 0.2100... 
                  0.2150 0.2200 0.2250 0.2300 0.2350 0.2400...
                  0.2450]...
                  + r_shift;


n_z = length(z_probe_pcb);

%　列方向が同じr位置にそろえる
% ch1-25, ch26-50, ch51-76, ch77-101, ch102-126 ch63は飛ばすから注意

%% bz signal former f1
f1=figure(Visible=figure_switch(1));
f1.WindowState = 'maximized';
for i=1:r
    for j=1:col1
            subplot(r,col1,(i-1)*col1+j)
            if ok_bz(r_ch*(i-1)+j)==1 %okなチャンネルはそのままプロット
                plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j))
                % hold on
                % plot(t_start:t_end,bz_s(t_start:t_end,r_ch*(i-1)+j))
            else %NGなチャンネルは赤色点線でプロット
                plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j),'r:')
            end
            title(num2str(p_ch(i,j)));
            %xlim([t_start t_end]);
            xticks([t_start t_end]);
            ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bz signal, r = 0.06-0.18')

% save plot image
foldername = strcat(pathname.fig, num2str(date));
if exist(foldername, 'dir') == 0
    mkdir(foldername);
end
savename = strcat(foldername, '/shot', num2str(shot(1), '%04i'), '_', 'bz_signal1.png');
exportgraphics(gcf,savename, 'Resolution',300);


%% bz signal latter f2
f2=figure(Visible=figure_switch(2));
f2.WindowState = 'maximized';
for i=1:r
    for j=col1+1:col1+col2
            subplot(r,col2,(i-1)*col2+j-col1)
            if ok_bz(r_ch*(i-1)+j)==1
                plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j))
                % hold on
                % plot(t_start:t_end,bz_s(t_start:t_end,r_ch*(i-1)+j))
            else
                plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j),'r:')
            end
            title(num2str(p_ch(i,j)));
            %xlim([t_start t_end]);
            xticks([t_start t_end]);
            ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bz signal, r = 0.185-0.245')

% save plot image
savename = strcat(foldername, '/shot', num2str(shot(1), '%04i'), '_', 'bz_signal2.png');
exportgraphics(gcf,savename, 'Resolution',300);


%% 横軸z, 縦軸Bzのプロット
f3=figure(Visible=figure_switch(3));
f3.WindowState = 'maximized';
styles = ["-*r","-^r","-or","-squarer","-diamondr",...
    "-*b","-^b","-ob","-squareb","-diamondb",...
    "-*k","-^k","-ok","-squarek","-diamondk",...
      "-*c","-^c","-oc","-squarec","-diamondc",...
        "-*m","-^m","-om","-squarem","-diamondm"];
tiles = tiledlayout(2,1);
sgtitle(strcat('t=',num2str(t),' us'))

nexttile
hold on


for i=1:25
    % rpos = 0.06[m]の位置のbzの値 zline == [1,26,51,77,102]
    zline=(1:25:n_z*25-24)+(i-1);
    %zline = [1,26,51,77,102];
    % r = 0.21mの信号
    %zline = [18,43,69,94,119];
    bz_zline=bz(t,zline);
    bz_zline(ok_bz(zline)==false)=NaN;
    plot(z_probe_pcb,bz_zline,styles(i),'MarkerSize',12)
    clear bz_zline
end
hold off
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10', ...
    'r11','r12','r13','r14','r15','r16','r17','r18','r19','r20',...
    'r21','r22','r23','r24','r25' ...
    ,Location='eastoutside')
title('Before interpolation')

Bz_interped = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
nexttile
hold on
for i=1:25
    % rpos = 0.06[m]の位置のbzの値 zline == [1,26,51,77,102]
    zline=(1:25:n_z*25-24)+(i-1);
    %zline = [1,26,51,77,102];
    % r = 0.21mの信号
    %zline = [18,43,69,94,119];
    bz_zline=Bz_interped(zline);
    plot(linspace(min(zpos_bz),max(zpos_bz),n_z),bz_zline,styles(i),'MarkerSize',12)
    clear bz_zline
end
hold off
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
%ylim([0 0.12])
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10', ...
    'r11','r12','r13','r14','r15','r16','r17','r18','r19','r20',...
    'r21','r22','r23','r24','r25' ...
    ,Location='eastoutside')
title('After interpolation')
tiles.TileSpacing = 'compact';
tiles.Padding = 'compact';


%% 横軸r,縦軸Bzのプロット(z = 0, -0.021)(x点は負の方向にずれていることが多い)
f4=figure(Visible=figure_switch(4));
f4.WindowState = 'maximized';

tile = tiledlayout(3,2);

for probe_num = 1:r
    rline = [(1+(probe_num-1)*r_ch):(probe_num*r_ch)]; 
    bz_rline=bz(t,rline);
    bz_rline(ok_bz(rline)==false)=NaN;

    nexttile
    plot(r_probe_pcb,bz_rline,'r.','MarkerSize',12)
    clear bz_rline r_line
    title(strcat('z = ',num2str(z_probe_pcb(probe_num))));
    yline(0,'k--')
    %xline(r_probe_pcb,':','LineWidth',1.5)
    ylim([-0.02,0.1])
    xlim([min(r_probe_pcb),max(r_probe_pcb)])
end
tiles.TileSpacing = 'compact';
tiles.Padding = 'compact';
title(tile,strcat('t=',num2str(t),' us'))
xlabel(tile,'r[m]')
ylabel(tile,'bz[T]')
% save plot image
savename = strcat(foldername, '/shot', num2str(shot(1), '%04i'), '_', 'bz_z.png');
exportgraphics(gcf,savename, 'Resolution',300);


%%
hidden = find(figure_switch == 'off');
figures = [f1,f2,f3,f4];
for i = hidden
    close(figures(i));
end
%}
