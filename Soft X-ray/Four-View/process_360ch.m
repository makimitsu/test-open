function [data2D,grid2D]=process_360ch(date, IDXlist, tfshot, pathname, n,i_EF,trange)
%%%%%%%%%%%%%%%%%%%%
% ポロイダル磁気面の時間発展プロット
%%%%%%%%%%%%%%%%%%%%
pathname.processed_data = "G:\My Drive\X-ray\Local Code Files";%processed data の保管場所
pathname.fig_psi = "G:\My Drive\X-ray\Local Code Files"; % plot画像の保存先
pathname.fig_psi = "G:\My Drive\X-ray\Data\PSIOUT"; % plot画像の保存先

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,date);
n_data=numel(IDXlist);% 計測データ数
shotlist_a039 = T.a039(IDXlist);
shotlist = shotlist_a039;
colorplot = 'Jt'; % 【input】color plot parameter
fitrate_plot_mode = 1; % 【input】fitrate plot parameter
sheet_plot_mode = 0;% 【input】fitrate plot parameter
start = 462;
dt = 1;

for shot = shotlist
filename = strcat(pathname.processed_data,'\a03839_',num2str(shot),'.mat');

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
% 図面の作成
% figureウインドウを画面左下隅から右に$1ピクセル、上に$2ピクセルの位置に配置 幅$3ピクセル、高さ$4ピクセル
%figure('Position',[0 0 1600 900],'visible','on');
figure('Position',[0 0 800 600],'visible','on');

total = 6;

tile = tiledlayout(2,3); % Figureを4行4列に分割。左上から右向きm番目の位置に図を描画
for m=1:total
    i = start + (m-1).*dt + 1 - 400;
    t = data2D.trange(i); % trange = 400:600
    nexttile
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),-data2D.Jt(:,:,i),30,'LineStyle','none'); clim([-0.7e+6, 0.7*1e+6]); % Jt
    label = 'Jt [A/m^2]';
    hold on
    % Jt等高線の表示
    contour(grid2D.zq(1,:),grid2D.rq(:,1),-data2D.Jt(:,:,i),[2e5:1e5:10e5],'white','LineWidth',0.9)
    hold off

    colormap(jet)
    axis image % 各軸のデータ単位を合わせる
    axis tight manual % 軸の範囲をデータ範囲と合わせて固定

    hold on
    % psi等高線の表示
    contour(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),-20e-3:0.3e-3:40e-3,'black','LineWidth',0.5) % psi
    xlim([-0.0525 0.0525])    
    xticks([-0.05 0 0.05])
    xticklabels({'-0.05', '0', '0.05'})
    xtickangle(0)
    xlabel('z[m]','FontSize',10)
    ylabel('r[m]','FontSize',10)
    hold off
    title(string(t) + 'us','FontSize',10)
end
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = label;
tile.TileSpacing = 'compact'; % 各tileの間隔を縮める
tile.Padding = 'compact'; % 各title内の余白を縮める

fontsize(gcf,25,'pixels')
end