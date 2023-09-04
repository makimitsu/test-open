%%%%%%%%%%%%%%%%%%%%%%%%
%200ch用新規pcbプローブと装置の磁場信号極性チェック
%dtacqのshot番号を直接指定する場合
%%%%%%%%%%%%%%%%%%%%%%%%

%figure_switch = ["on","on","off","off","on","on"];%bz1,bz2,bt1,bt2,bz_vs_z,bt_vs_z

r = 7;%プローブ本数＝グラフ出力時の縦に並べる個数
col = 10;%グラフ出力時の横に並べる個数
y_upper_lim = 0.03;%縦軸プロット領域（b_z上限）
y_lower_lim = -0.03;%縦軸プロット領域（b_z下限）
t_start=1;%横軸プロット領域（開始時間）
t_end=1000;%横軸プロット領域（終了時間）
t = 470;
r_shift = 0;
% z_probe_pcb = [-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17];
z_probe_pcb = [-0.2975 -0.255 -0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17 0.255 0.2975];
r_probe_pcb    = [0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33]+r_shift;
n_z = length(z_probe_pcb);
% r_ch=col1+col2;%r方向から挿入した各プローブのチャンネル数

f6=figure(Visible=figure_switch(6));
f6.WindowState = 'maximized';
%styles = ["-*","-^","-v","-<","->","-o","-square","-diamond","-pentagram","-hexagram"];
styles
tiles = tiledlayout(2,1);
sgtitle(strcat('t=',num2str(t),' us'))
nexttile
hold on
%for i=1:n_z
    %if i == 6
    i = 6;
    rline= i*10-9:i*10;%14本のプローブの一番内側のbz点：1,11,21,31,41,51,61,71,81,91,101,111,121,131
    % Xpointが見えそうな中心付近のプローブ（-0.0315 -0.0105 0.0105 0.0315）のr-bz plotを書く
    % z=-0.0315: bz(:,51:60)→plot
    % z=-0.0105: bz(:,61:70)→plot
    % z= 0.0105: bz(:,71:80)→plot
    % z= 0.0315: bz(:,81:90)→plot

    bz_rline=bz(t,rline);
    bz_rline(ok_bz(rline)==false)=NaN;
    plot(r_probe_pcb,bz_rline,'r.','MarkerSize',12)
    clear bz_zline
    %end
    hold off
    yline(0,'k--')
    xline(z_probe_pcb,':','LineWidth',1.5)
    %legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
    legend(strcat(z_probe_pcb(i)))
    title('Before interpolation')
%end
% hold off
% yline(0,'k--')
% xline(z_probe_pcb,':','LineWidth',1.5)
% %legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
% legend(strcat(z_probe_pcb(i)))
% title('Before interpolation')

%{
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
%}
tiles.TileSpacing = 'compact';
tiles.Padding = 'compact';
