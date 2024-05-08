%%% cal & plot ExB drift velocity %%%
function ESPdata2D = cal_ESP(pathname,ESP)
n_ch = 21;%静電プローブCH数
res_ratio = 50;%静電プローブ分圧比

savename = [pathname.mat,'/ESP/mesh',num2str(ESP.mesh),'_',num2str(ESP.date),'_shot',num2str(ESP.shotlist(1)),'-',num2str(ESP.shotlist(end)),'.mat'];
if exist(savename,"file")
    load(savename,'ESPdata2D')
else
    z = linspace(-0.15,0.15,ESP.mesh);%プロットメッシュZ座標[m]
    z_probe = linspace(-0.15,0.15,21);%静電プローブ計測点Z座標[m]
    ESPdata2D.zprobe = z_probe;
    ng_ch = [4 6 16];%死んだCH
    z_probe(ng_ch) = [];
    r = linspace(min(ESP.rlist),max(ESP.rlist),ESP.mesh)*1E-3;%プロットメッシュR座標[m]
    r_probe = unique(ESP.rlist)*1E-3;%静電プローブ計測点R座標[m]
    ESPdata2D.rprobe = r_probe;

    phi = zeros(numel(ESP.trange),n_ch,numel(r_probe));
    cnt_r = zeros(numel(r_probe),1);
    for i = 1:numel(ESP.shotlist)
        idx_r = find(r_probe==ESP.rlist(i)*1E-3);
        cnt_r(idx_r) = cnt_r(idx_r) + 1;
        filename = sprintf("%s%03d%s",[pathname.ESP '/' num2str(ESP.date) '/ES_' num2str(ESP.date)], ESP.shotlist(i), '.csv');
        ESPdata = readmatrix(filename,'Range',sprintf('B%d:V%d',ESP.trange(1)*10+2,ESP.trange(end)*10+2));
        phi(:,:,idx_r) = (phi(:,:,idx_r)*(cnt_r(idx_r)-1) + ESPdata)/cnt_r(idx_r);
    end
    phi(:,ng_ch,:) = [];%死んだCHを除去
    phi = fliplr(phi);%(CH1のZ座標)>(CH2のZ座標)>...のため、列を反転
    phi = phi.*res_ratio;%分圧比を掛ける

    [ESPdata2D.zq,ESPdata2D.rq] = meshgrid(z,r);
    ESPdata2D.trange = ESP.trange;
    ESPdata2D.phi = zeros(numel(ESP.trange),ESP.mesh,ESP.mesh);
    ESPdata2D.Ez = zeros(numel(ESP.trange),ESP.mesh,ESP.mesh);
    ESPdata2D.Er = zeros(numel(ESP.trange),ESP.mesh,ESP.mesh);
    for i = 1:numel(ESP.trange)
        ESPdata2D.phi(i,:,:) = griddata(z_probe,r_probe,squeeze(phi(i,:,:))',ESPdata2D.zq,ESPdata2D.rq);
        ESPdata2D.Ez(i,:,2:end) = -diff(squeeze(ESPdata2D.phi(i,:,:)),1,2)/(z(2)-z(1));
        ESPdata2D.Er(i,2:end,:) = -diff(squeeze(ESPdata2D.phi(i,:,:)),1,1)/(r(2)-r(1));
    end
    save(savename,'ESPdata2D')
end
