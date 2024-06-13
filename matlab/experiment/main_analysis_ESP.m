%%% 静電プローブデータ %%%
% close all
clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common';
run define_path.m

%-----------【input】-------------
ESP.date = 230830;%【input】計測日
ESP.shotlist = [11 13 14 16 17 20 22:26 28 29 32 34 37 41:45 47 49:51 54:57 59 60];%【input】静電プローブ解析shotlist(同一オペレーション)
ESP.trange = 450:0.1:500;%【input】計算時間範囲(0.1刻み)
plot_time = 470;%【input】プロット時刻[us]
analysis_axis = 'phi';%【input】解析横軸'Z','R','ZR'
ESP.ng_ch = [4 6 16 21];%【input】死んだCH番号
ESP.ng_shot = [14 17 24 25 28 43 44 51 55 57 60 50 56];%【input】除くshot番号14 17 24 25 28 43 44 51 55 57 60

%------定数------
ESP.n_ch = 21;%静電プローブCH数
ESP.ch = linspace(1,ESP.n_ch,ESP.n_ch);%静電プローブCH番号
ESP.res_ratio = 50;%静電プローブ分圧比
ESP.z_probe = linspace(0.15,-0.15,21);%%静電プローブ計測点Z座標[m](CH1がZ=0.15m、CH21が-0.15mであることに注意!)

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,ESP.date);
ESP.rlist=T.ESProbeRPosition_mm_(ESP.shotlist);%静電プローブr座標[mm]
ESP.rlist=ESP.rlist*1E-3;%[mm]を[m]に変換

%ファイル名をshot番号リストに対応して命名
savename = [pathname.mat,'/ESP_analysis/',num2str(ESP.date),'_shot'];
for i_shot = 1:numel(ESP.shotlist)
    if i_shot == 1
        savename =[savename,num2str(ESP.shotlist(i_shot))];
    else
        if ESP.shotlist(i_shot) == ESP.shotlist(i_shot-1)+1%連番の場合間の番号をファイル名に含まない
            if i_shot < numel(ESP.shotlist)
                if ESP.shotlist(i_shot+1) > ESP.shotlist(i_shot)+1
                    savename =[savename,'-',num2str(ESP.shotlist(i_shot))];
                end
            else
                savename =[savename,'-',num2str(ESP.shotlist(i_shot))];
            end
        else%連番でない場合ファイル名に含む
            savename =[savename,'_',num2str(ESP.shotlist(i_shot))];
        end
    end
end
savename = [savename,'.mat'];

if exist(savename,"file")
    load(savename,'ESPdata2D')
else
    ESPdata2D.trange = ESP.trange;
    ESPdata2D.phi = zeros(numel(ESP.trange),ESP.n_ch,numel(ESP.shotlist));%生信号
    ESPdata2D.LPF_phi = zeros(size(ESPdata2D.phi));%Low pass filtered
    for i_shot = 1:numel(ESP.shotlist)
        filename = sprintf("%s%03d%s",[pathname.ESP '/' num2str(ESP.date) '/ES_' num2str(ESP.date)], ESP.shotlist(i_shot), '.csv');
        ESPdata = readmatrix(filename,'Range',sprintf('B%d:V%d',ESP.trange(1)*10+2,ESP.trange(end)*10+2));
        ESPdata2D.phi(:,:,i_shot) = ESPdata;
        ESPdata2D.phi = ESPdata2D.phi.*ESP.res_ratio;%分圧比を掛ける
        for i_ch = 1:size(ESPdata2D.phi,2)
            ESPdata2D.LPF_phi(:,i_ch,i_shot) = lowpass(ESPdata2D.phi(:,i_ch,i_shot),1E-4);
        end
    end
    save(savename,'ESPdata2D')
end

%死んだCHを除去
if ESP.ng_ch
    ESP.ch(ESP.ng_ch) = [];
    ESP.z_probe(ESP.ng_ch) = [];
    ESPdata2D.phi(:,ESP.ng_ch,:) = [];
    ESPdata2D.LPF_phi(:,ESP.ng_ch,:) = [];
end
%除きたいshotを除去
if ESP.ng_shot
    idx_ng_shot = knnsearch(ESP.shotlist',ESP.ng_shot');
    ESP.shotlist(idx_ng_shot) = [];
    ESP.rlist(idx_ng_shot) = [];
    ESPdata2D.phi(:,:,idx_ng_shot) = [];
    ESPdata2D.LPF_phi(:,:,idx_ng_shot) = [];
end
idx_time = knnsearch(ESPdata2D.trange',plot_time);
phi_at_t = transpose(squeeze(ESPdata2D.LPF_phi(idx_time,:,:)));
switch analysis_axis
    case 'Z'
        figure('Position', [0 0 1500 1500],'visible','on');
        zz = linspace(min(ESP.z_probe),max(ESP.z_probe),100);
        mov_z_phi_at_t =  movmean(phi_at_t,5,2);
        for i_shot = 1:numel(ESP.shotlist)
            yy = spline(ESP.z_probe,mov_z_phi_at_t(i_shot,:),zz);
            subplot(6,6,i_shot)
            plot(ESP.z_probe,phi_at_t(i_shot,:),'bo',zz,yy,'r-')
            title(['Shot',num2str(ESP.shotlist(i_shot)),': R=',num2str(ESP.rlist(i_shot)),'mm'])
            xlabel('Z [m]')
            ylabel('Floating Potential [V]')
            ylim([-200 250])
        end
        sgtitle([num2str(plot_time),'us'])
    case 'R'
        figure('Position', [0 0 1500 1500],'visible','on');
        rr = linspace(min(ESP.rlist),max(ESP.rlist),100);
        mov_r_phi_at_t =  movmean(phi_at_t,5,1);
        for i_ch = 1:numel(ESP.ch)
            yy = spline(ESP.rlist,mov_r_phi_at_t(:,i_ch),rr);
            subplot(4,5,i_ch)
            plot(ESP.rlist,phi_at_t(:,i_ch),'bo',rr,yy,'r-')
            title(['CH',num2str(ESP.ch(i_ch)),': Z=',num2str(ESP.z_probe(i_ch)),'m'])
            xlabel('R [m]')
            ylabel('Floating Potential [V]')
            ylim([-200 250])
        end
        sgtitle([num2str(plot_time),'us'])
    case 'ZR'
        figure('Position', [0 0 1500 1500],'visible','on');
        [mesh_z,mesh_r] = meshgrid(ESP.z_probe,ESP.rlist);
        zz = linspace(min(ESP.z_probe),max(ESP.z_probe),100);
        rr = linspace(min(ESP.rlist),max(ESP.rlist),100);
        mov_zr_phi_at_t =  movmean(phi_at_t,3,1);
        mov_zr_phi_at_t =  movmean(mov_zr_phi_at_t,3,2);
        grid_phi = griddata(ESP.z_probe,ESP.rlist',mov_zr_phi_at_t,zz',rr);
        surf(zz,rr,grid_phi)
        hold on
        stem3(mesh_z,mesh_r,phi_at_t,'--*m')
        xlabel('Z [m]')
        ylabel('R [m]')
        c = colorbar;
        c.Label.String = 'Floating Potential [V]';
        clim([-200 250])
    case 'phi'
        figure('Position', [0 0 1500 1500],'visible','on');
        [mesh_z,mesh_r] = meshgrid(ESP.z_probe,ESP.rlist);
        mov_zr_phi_at_t =  movmean(phi_at_t,3,1);
        mov_zr_phi_at_t =  movmean(mov_zr_phi_at_t,3,2);
        contourf(mesh_z,mesh_r,mov_zr_phi_at_t,100,'edgecolor','none');
        c = colorbar;
        clim([-240 240])
        c.Label.String = 'Floating phi [V]';
        colormap(redblue(300))
        daspect([1 1 1])
        view([90 -90])%RZ反転
end