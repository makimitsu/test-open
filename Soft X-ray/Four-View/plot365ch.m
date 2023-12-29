%%%%%%%%%%%%%%%%%%%%
% ポロイダル磁気面の時間発展プロット
%%%%%%%%%%%%%%%%%%%%
clear all
pathname.processed_data = "G:\My Drive\X-ray\Data\PSIOUT\MATDATA";%processed data の保管場所
pathname.fig_psi = "G:\My Drive\X-ray\Data\PSIOUT\FIG"; % plot画像の保存先

shotlist = 3456;% FC合体shot7-13:3389-3395; SEP合体shot15-26:3450-3461
colorplot = 'Jt'; % 【input】color plot parameter
fitrate_plot_mode = 1; % 【input】fitrate plot parameter
sheet_plot_mode = 0;% 【input】fitrate plot parameter
start = 460;
dt = 1;

for shot = shotlist
%filename = strcat(pathname.processed_data,'\a039_',num2str(shot),'.mat');
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
figure('Position',[100 100 640 480],'visible','on');
set(gcf,'Name','shot' + string(shot(2)),'NumberTitle','off');

total = 6;

if sheet_plot_mode
    for m = 1:total
        grid2D = calc_sheet_width(grid2D, data2D,start,dt,total,m);
    end
end

tile = tiledlayout(2,3); % Figureを4行4列に分割。左上から右向きm番目の位置に図を描画
for m=1:total
    i = start + (m-1).*dt + 1 - 400;
    t = data2D.trange(i); % trange = 400:600
    nexttile



    
    % flag_Ip_area = 0;
    switch colorplot
        case 'Psi'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none'); clim([-8e-3,8e-3]); % psi
            label = 'Psi [Wb]';
        case 'Bz'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none'); clim([-0.1, 0.1]); % Bz
            label = 'Bz [T]';
        case 'Bt'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-50e-3:0.2e-3:50e-3,'LineStyle','none'); clim([-0.06, 0.06]); % Bt
            % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),0.05:0.01:0.3,'LineStyle','none'); clim([0.05, 0.3]); % TF only の場合
            label = 'Bt [T]';
        case 'Bt_ex'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt_ex(:,:,i),0.05:0.01:0.3,'LineStyle','none'); clim([0.05, 0.3]); % Bt_ex
            label = 'Bt_{ex} [T]';
        case 'Babs'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Babs(:,:,i),0:0.01:0.3,'LineStyle', 'none'); clim([0, 0.3]); % Babs
            label = '|B| [T]';
        case 'mg_P'
                contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.mg_P(:,:,i),1e2:0.1e2:1.5e4,'LineStyle', 'none'); clim([1e2, 1.5e4]); % Babs
                label = 'B^2/(2mu)';
        case 'Jt'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),-data2D.Jt(:,:,i),30,'LineStyle','none'); clim([-0.7e+6, 0.7*1e+6]); % Jt
            label = 'Jt [A/m^2]';
            hold on
            % Jt等高線の表示
            
            contour(grid2D.zq(1,:),grid2D.rq(:,1),-data2D.Jt(:,:,i),[2e5:2e5:10e5],'white','LineWidth',1)
            hold off
        case 'Jz'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Jz(:,:,i),30,'LineStyle','none'); clim([-0.8*1e+6, 0.8*1e+6]); % Jt
            label = 'Jz [A/m^2]';
             % Jt等高線の表示
            contour(grid2D.zq(1,:),grid2D.rq(:,1),-data2D.Jt(:,:,i),3.0e5:1.5e5:6e5,'white','LineWidth',1)
        case 'Jr'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Jr(:,:,i),40,'LineStyle','none'); %clim([-8e-3,8e-3]); 
            label = 'Psi [Wb]';
        case 'Fabs'
            %contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Fabs(:,:,i),40,'LineStyle','none'); %clim([-8e-3,8e-3]); 
            label = '|j×B| [N]';
            hold on
            q = quiver(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.Fz(:,:,i)),squeeze(data2D.Fr(:,:,i)),1.5);
            q.Color = "k";
            q.LineWidth = 2;
            hold off
        case 'Et'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none'); clim([-500, 400]); % Et
            label = 'Et [V/m]';
        case 'eta'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.eta(:,:,i),20,'LineStyle','none'); clim([-10, 200]); % eta
            label = 'eta[Ω m]';
        case 'q'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.q(:,:,i),10,'LineStyle','none'); clim([-2, 2]); % q
            label = 'q';
        case 'p'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.P(:,:,i),20,'LineStyle','none'); clim([0, 3000]); % P
            label = 'P [Pa]';
        case 'dpdr'
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.dPdr(:,:,i),10,'LineStyle','none'); %clim([-3000, 3000]); % dPdr
            label = 'dP/dr [Pa/m]';
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
                        hold on
                    end
                end
            end
            % hold off
            label = 'None';
        otherwise
            disp('error: incorrect input for colorplot')
            return
    end
    
    %colormap(whitejet)
    colormap(turbo)
    axis image % 各軸のデータ単位を合わせる
    axis tight manual % 軸の範囲をデータ範囲と合わせて固定

    hold on
    % psi等高線の表示
    contour(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),-20e-3:0.9e-3:40e-3,'black','LineWidth',1.5) % psi
    %yline(grid2D.xpos(2,i))

    if sheet_plot_mode == 1
        xline(grid2D.lower_half_width(i),'b','LineWidth',1)
        xline(grid2D.upper_half_width(i),'b','LineWidth',1)
        yline(grid2D.lower_half_height(i),'b','LineWidth',1)
        yline(grid2D.upper_half_height(i),'b','LineWidth',1)
        
    end
    
%    plot(grid2D.xpos(1,i),grid2D.xpos(2,i),'kx')%X点

    % 
    % 
    % for i=1:size(grid2D.zprobepcb_a39,2)
    %     for j=1:size(grid2D.rprobepcb_t,2)
    %         if grid2D.ok_bz_matrix_a39(j,i) == 1 % bz測定点の表示
    %             plot(grid2D.zprobepcb_a39(i),grid2D.rprobepcb_a39(j),'b.','MarkerSize',8);
    %         end
    %         if grid2D.ok_bt_matrix_a39(j,i) == 1 % bt測定点の表示
    %             plot(grid2D.zprobepcb(i),grid2D.rprobepcb_t(j),'k.','MarkerSize',6);
    %         end
    %     end
    % end
    % 
    % for i=1:size(grid2D.zprobepcb_a38,2)
    %     for j=1:size(grid2D.rprobepcb_a38,2)
    %         if grid2D.ok_bz_matrix_a38(j,i) == 1 % bz測定点の表示
    %             plot(grid2D.zprobepcb_a38(i),grid2D.rprobepcb_a38(j),'g.','MarkerSize',6);
    %         end
    %         if grid2D.ok_bt_matrix_a38(j,i) == 1 % bt測定点の表示
    %             plot(grid2D.zprobepcb(i),grid2D.rprobepcb_t(j),'k.','MarkerSize',6);
    %         end
    %     end
    % end

  
    % viscircles([0.2135,0.31],0.02,'Color','k');
    % viscircles([-0.2135,0.31],0.02,'Color','k');
    % viscircles([0.35,0.22],0.0375,'Color','k');
    % viscircles([-0.35,0.22],0.0375,'Color','k');
    xlim([-0.0525 0.0525])
    
    xticks([-0.05 0 0.05])
    % xticklabels({'-0.05', '0', '0.05'})
    xticklabels({})
    yticklabels({})
    xtickangle(0)
    daspect([0.8 1 1]);
    

    %xlim([-0.055 0.055])
    %ylim([0.06 0.245])
    %xlim([-0.13 0.13])
    

    % xlabel('z[m]','FontSize',10)
    % ylabel('r[m]','FontSize',10)
    hold off

    title(string(t) + 'us','FontSize',10)
end
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = label;
tile.TileSpacing = 'compact'; % 各tileの間隔を縮める
tile.Padding = 'compact'; % 各title内の余白を縮める

fontsize(gcf,25,'pixels')
%%% save plot image

foldername_psi = strcat(pathname.fig_psi, '/',num2str(date));
if exist(foldername_psi, 'dir') == 0
   mkdir(foldername_psi);
end
savename_psi = strcat(foldername_psi, '/shot', num2str(shot(2), '%04i'), '_', colorplot, '_unified_with_axis_',num2str(start),'to',num2str(t),'.png');
savename_psi_fig = strcat(foldername_psi, '/shot', num2str(shot(2), '%04i'), '_', colorplot, '_unified_with_axis_',num2str(start),'to',num2str(t),'.fig');
exportgraphics(gcf,savename_psi, 'Resolution',300);
saveas(gcf,savename_psi_fig);
%% plot sheet width time change





end