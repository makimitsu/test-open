%%% cal & plot ExB drift velocity %%%
function data2D = cal_ESP(pathname,date,IDXlist,rlist,trange,n_mesh)
n_ch = 21;%静電プローブCH数
res_ratio = 50;%静電プローブ分圧比

z = linspace(-0.15,0.15,n_mesh);%プロットメッシュZ座標[m]
z_probe = linspace(-0.15,0.15,21);%静電プローブ計測点Z座標[m]
ng_ch = [4 6 16];%死んだCH
z_probe(ng_ch) = [];
r = linspace(min(rlist),max(rlist),n_mesh)*1E-3;%プロットメッシュR座標[m]
r_probe = unique(rlist)*1E-3;%静電プローブ計測点R座標[m]

phi = zeros(numel(trange),n_ch,numel(r_probe));
cnt_r = zeros(numel(r_probe),1);
for i = 1:numel(IDXlist)
    idx_r = find(r_probe==rlist(i)*1E-3);
    cnt_r(idx_r) = cnt_r(idx_r) + 1;
    filename = sprintf("%s%03d%s",[pathname.ESP '/' num2str(date) '/ES_' num2str(date)], IDXlist(i), '.csv');
    ESPdata = readmatrix(filename,'Range',sprintf('B%d:V%d',trange(1)*10+2,trange(end)*10+2));
    phi(:,:,idx_r) = (phi(:,:,idx_r)*(cnt_r(idx_r)-1) + ESPdata)/cnt_r(idx_r);
end
phi(:,ng_ch,:) = [];%死んだCHを除去
phi = fliplr(phi);%(CH1のZ座標)>(CH2のZ座標)>...のため、列を反転
phi = phi.*res_ratio;%分圧比を掛ける

[data2D.phi_mesh_z,data2D.phi_mesh_r] = meshgrid(z,r);
data2D.phi_grid = zeros(numel(trange),n_mesh,n_mesh);
data2D.Ez_grid = zeros(numel(trange),n_mesh,n_mesh);
data2D.Er_grid = zeros(numel(trange),n_mesh,n_mesh);
for i = 1:numel(trange)
    data2D.phi_grid(i,:,:) = griddata(z_probe,r_probe,squeeze(phi(i,:,:))',data2D.phi_mesh_z,data2D.phi_mesh_r);
    data2D.Ez_grid(i,:,2:end) = -diff(squeeze(data2D.phi_grid(i,:,:)),1,2)/(z(2)-z(1));
    data2D.Er_grid(i,2:end,:) = -diff(squeeze(data2D.phi_grid(i,:,:)),1,1)/(r(2)-r(1));
end