NofCH = 24;%チャンネル数
NT = 9;%計測トリガ時刻数
NData = 5;%同一計測トリガでの計測ショット数
nz = 1;%z方向計測点数
factor = 0.3;%フロー矢印の長さ倍率
r_measured = zeros(NofCH/4,nz);%ベクトルプロットr座標1列目data1、2列目data2
z_measured = zeros(NofCH/4,nz);%ベクトルプロットz座標1列目data1、2列目data2
z_measured(:,1) = 2.1;%data1のz座標
r_measured(:,1) = [22.5, 25, 27.5, 30, 32.5, 35];%data1のr座標

fig = figure;
frames(NT) = struct('cdata', [], 'colormap', []); % 各フレームの画像データを格納する配列

for t = 1:NT
    Ti_multi = zeros(NofCH/4,NData);
    Vz_multi = zeros(NofCH/4,NData);
    Vr_multi = zeros(NofCH/4,NData);
    for ndata = 1:NData
        load(['magprobe/mat/save_',num2str(463+2*t),'us_',num2str(ndata),'.mat'],'V','Ti')
        Ti_multi(:,ndata) = Ti(1:NofCH/4,1);
        Vz_multi(:,ndata) = V(1:NofCH/4,1);
        Vr_multi(:,ndata) = V(1:NofCH/4,2);
    end
    Ti_mean = trimmean(Ti_multi,60,2);
    Vz_mean = trimmean(Vz_multi,60,2);
    Vr_mean = trimmean(Vr_multi,60,2);
    absV_mean = zeros(NofCH/4,1);
    if t == 1
        Ti_offset = Ti_mean;
    end
    for i = 1:NofCH/4
        absV_mean(i,1) = sqrt(Vz_mean(i,1)^2 + Vr_mean(i,1)^2);
    end
    
    %平均流速、平均温度をプロット
    %     Ticon = repmat(Ti_mean,1,2);%等高線図用平均イオン温度
    Ticon = repmat(Ti_mean - Ti_offset,1,2);%等高線図用平均イオン温度上昇
    zcon = z_measured(1,1);
    s = pcolor([zcon-0.8 zcon+0.8],r_measured,Ticon);
    s.FaceColor = 'interp';
    s.FaceAlpha = 0.6;%カラーマップ不透明度(0~1)
    s.EdgeAlpha = 0;
    colormap('jet')
    %     caxis([55 110])
    caxis([0 50])
    hold on
    plot(z_measured,r_measured,'xr','MarkerSize',5,'LineWidth',2);
    hold on
    if nz == 1
        q = quiver(z_measured,r_measured,Vz_mean*factor,Vr_mean*factor);
    end
    if nz == 2
        q = quiver(z_measured,r_measured,[V(:,1),V(:,3)]*factor,[V(:,2),V(:,4)]*factor);
    end
    q.LineWidth = 3;
    q.MaxHeadSize = 15;
    q.AutoScale = 'off';
    q.Color = 'black';
    xlim([min(z_measured,[],'all')-2.5 max(z_measured,[],'all')+2.5])
    ylim([min(r_measured,[],'all')-1 max(r_measured,[],'all')+1])
    %     title({[num2str(463+2*t),'us'];'Ion Flow [km/s]'},'Color','black','FontWeight','bold')
    title([num2str(463+2*t),'us'],'Color','black','FontWeight','bold')
    absV_mean = round(absV_mean,1);
    for j = 1:nz
        for i = 1:NofCH/4
            txt = text(z_measured(i,j)+0.9,r_measured(i,j)-0.2,num2str(absV_mean(i,j)));
            txt.FontSize = 16;
            txt.Color = 'k';
            txt.FontWeight = 'bold';
        end
    end
    xlabel('Z [cm]')
    ylabel('R [cm]')
    ax = gca;
    ax.FontSize = 12;
    grid on
    daspect([1 1 1])
    c = colorbar;
    % c.Label.String = 'Ion Temperature [eV]';
    c.Label.String = 'Increase of Ion Temperature [eV]';
    c.FontSize = 12;
    drawnow; % 描画を確実に実行させる
    frames(t) = getframe(fig); % 図を画像データとして得る
    hold off
end

%gifファイルを保存
filename = 'movie_flow.gif'; % ファイル名
for t = 1:NT
    [A, map] = rgb2ind(frame2im(frames(t)), 256); % 画像形式変換
    if t == 1
        imwrite(A, map, filename, 'gif', 'DelayTime', 1); % 出力形式(FPS)を設定
    else
        imwrite(A, map, filename, 'gif', 'DelayTime', 1, 'WriteMode', 'append'); % 2フレーム目以降は"追記"の設定も必要
    end
end
