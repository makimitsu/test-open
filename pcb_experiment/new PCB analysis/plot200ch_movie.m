%%%%%%%%%%%%%%%%%%%%
% ポロイダル磁気面のアニメーション作成
% （参考）https://qiita.com/tose2125/items/03d9ce40b1b6bde7d36d
%%%%%%%%%%%%%%%%%%%%
clear all
pathname.processed_data="C:\Users\uswk0\OneDrive\ドキュメント\data\pre_processed";%processed data の保管場所
pathname.fig_psi = "C:\Users\uswk0\OneDrive\ドキュメント\data\figure"; % plot画像の保存先

shot = 2113;        % 【input】shot number of a039
colorplot = 'Psi';  % 【input】color plot parameter
exporttype = 'mp4'; % 【input】'mp4' or 'gif'
t_start = 460;      % 【input】start time
t_end = 650;        % 【input】end time

filename = strcat(pathname.processed_data,'/a039_',num2str(shot),'.mat');


if exist(filename,'file') == 0
    disp(['File:',filename,' does not exist']);
    return
end
load(filename);

if isstruct(grid2D) == 0 || isstruct(data2D) == 0
    disp('data incomplete/corrupted');
    return
end

%%%%%%%%%%%%%%%%%%%%%
% アニメーションの作成

fig = figure; % Figure オブジェクトの生成
fig.Color = 'white'; % 背景色を白に設定
fig.ToolBar = 'none'; % ツールバーを非表示

len = t_end - t_start + 1;
idx_start = t_start - data2D.trange(1) + 1;
idx_end = idx_start + len - 1;
if idx_start < 1 || idx_end > size(data2D.trange,2) || idx_start > idx_end
    disp('error: incorrect input for trange')
end
frames(len) = struct('cdata', [], 'colormap', []); % 各フレームの画像データを格納する配列

idx_flame = 1;
for i = idx_start:idx_end
    % % psi等高線の表示
    % contour(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),-20e-3:0.3e-3:40e-3,'black','LineWidth',1) % psi
    % hold on
    switch colorplot
        case 'Psi'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none'); clim([-6e-3,6e-3]); % psi
            label = 'Psi [Wb]';
        case 'Bz'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none'); clim([-0.1, 0.1]); % Bz
            label = 'Bz [T]';
        case 'Bt'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-50e-3:0.2e-3:50e-3,'LineStyle','none'); clim([-0.08, 0.08]); % Bt
            label = 'Bt [T]';
        case 'Bt_ex'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt_ex(:,:,i),0:0.01:0.4,'LineStyle','none'); clim([0, 0.4]); % Bt_ex
            label = 'Bt_{ex} [T]';
        case 'Babs'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Babs(:,:,i),0:0.04:0.5, 'LineStyle', 'none'); clim([0, 0.5]); % Babs
            label = '|B| [T]';
        case 'Jt'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none'); clim([-0.8*1e+6, 0.8*1e+6]); % Jt
            label = 'Jt [A/m^2]';
        case 'Et'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none'); clim([-500, 400]); % Et
            label = 'Et [V/m]';
        case 'Ip_area'
            zpos = grid2D.zq(1,:);
            rpos = grid2D.rq(:,1);
            %%% セパラトリクスとなるpsiを求める
            max_psi = zeros(size(zpos));
            for j = 1:size(zpos,2)
                max_psi(j) = max(data2D.psi(:,j,i));
            end
            max_psi_zng = max_psi(1:size(zpos,2)/2);
            max_psi_zpos = max_psi(size(zpos,2)/2+1:size(zpos,2));
            min_psi_zng = min(max_psi_zng(max_psi_zng>=0));
            min_psi_zpos = min(max_psi_zpos(max_psi_zpos>=0));
            for j = 1:size(zpos,2)
                if max_psi(j) == min_psi_zng
                    z_sep_idx_ng = j;
                elseif max_psi(j) == min_psi_zpos
                    z_sep_idx_pos = j;
                end
            end
            psi_sep = max(min_psi_zng,min_psi_zpos);
            for k = 1:size(rpos,1)
                r = rpos(k);
                for j = 1:size(zpos,2)
                    z = zpos(j);
                    if data2D.psi(k,j,i) >= psi_sep && j>=z_sep_idx_ng && j <= z_sep_idx_pos
                        plot(z,r,'r.','MarkerSize',6);
                        % hold on
                    end
                end
            end
            hold off
            label = 'None';
        otherwise
            disp('error: incorrect input for colorplot')
            return
    end
    colormap(whitejet)
    axis image % 各軸のデータ単位を合わせる
    axis tight manual % 軸の範囲をデータ範囲と合わせて固定
    hold on
    % psi等高線の表示
    contour(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),-20e-3:0.3e-3:40e-3,'black','LineWidth',1) % psi
    hold off
    title(string(data2D.trange(i)) + 'us')
    cb = colorbar;
    cb.Label.String = label;

    drawnow; % 描画を確実に実行させる
    frames(idx_flame) = getframe(fig); % 図を画像データとして得る
    idx_flame = idx_flame + 1;
end


%%% save movie
foldername = strcat(pathname.fig_psi, num2str(date),'/',num2str(shot(1)));
if exist(foldername, 'dir') == 0
    mkdir(foldername);
end
savename = strcat(foldername, '/shot', num2str(shot(1), '%04i'), '_', colorplot);


switch exporttype
    case 'mp4'
        video = VideoWriter(strcat(savename, '.mp4'), 'MPEG-4'); % ファイル名と出力形式を設定
        open(video); % 書き込むファイルを開く
        writeVideo(video, frames); % ファイルに書き込む
        close(video); % 書き込むファイルを閉じる
    case 'gif'
        filename = strcat(savename, '.gif'); % ファイル名
        for i = 1:len
            [A, map] = rgb2ind(frame2im(frames(i)), 256); % 画像形式変換
            if i == 1
                imwrite(A, map, filename, 'gif', 'DelayTime', 1/30); % 出力形式(30FPS)を設定
            else
                imwrite(A, map, filename, 'gif', 'DelayTime', 1/30, 'WriteMode', 'append'); % 2フレーム目以降は"追記"の設定も必要
            end
        end
    otherwise
        disp('error: incorrect input for exportype')
        return
end

