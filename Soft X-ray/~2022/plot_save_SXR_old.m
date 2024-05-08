function plot_save_SXR_old(B_z,r_probe,z_probe,range,date,shot,t,EE1,EE2,show_localmax,show_xpoint,save,filter,NL)
% plot SXR emission and psi in rz plane and save the result
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   1d array of integer: range, [zmin1,zmax1,zmin2,zmax2,rmin,rmax]
%   integer: date, date of the experiment
%   integer: shot, shot number
%   integer: t, time of interest (us)
%   2d array of double: EE1,EE2, local distribution of SXR emission
%   boolean: show_xpoint, option for showing the x-point
%   boolean: show_localmax, option for showing the local maximum point
%   boolean: save, option for saving the reconstruction result
%   boolean: filter, option for applying non-linear mean (NLM) filter
%   boolean: NL, option for using non-linear reconstruction

% 表示範囲の設定に使うパラメータを取得
range = range./1000;
zmin1 = range(1);
zmax1 = range(2);
zmin2 = range(3);
zmax2 = range(4);
rmin = range(5);
rmax = range(6);
rmin2 = range(5);
rmax2 = range(6);
% rmin2 = range(5)-0.01;
% rmax2 = range(6)-0.01;
r_space_SXR1 = linspace(rmin,rmax,size(EE1,1));
r_space_SXR2 = linspace(rmin2,rmax2,size(EE2,1));
z_space_SXR1 = linspace(zmin1,zmax1,size(EE1,2));
z_space_SXR2 = linspace(zmin2,zmax2,size(EE2,2));

r_range1 = find(0.060<=r_space_SXR1 & r_space_SXR1<=0.300);
% r_range1 = find(0.060<=r_space_SXR1 & r_space_SXR1<=0.330);
r_space_SXR1 = r_space_SXR1(r_range1);
r_range2 = find(0.060<=r_space_SXR2 & r_space_SXR2<=0.300);
% r_range2 = find(0.060<=r_space_SXR2 & r_space_SXR2<=0.330);
r_space_SXR2 = r_space_SXR2(r_range2);
z_range1 = find(-0.20<=z_space_SXR1 & z_space_SXR1<=0.20);
z_range2 = find(-0.20<=z_space_SXR2 & z_space_SXR2<=0.20);

z_space_SXR1 = z_space_SXR1(z_range1);
z_space_SXR2 = z_space_SXR2(z_range2);

EE1 = EE1(r_range1,z_range1);
EE2 = EE2(r_range2,z_range2);

z_space = linspace(z_probe(1),z_probe(end),50);
r_space = linspace(r_probe(1),r_probe(end),50);
[psi_mesh_z,psi_mesh_r] = meshgrid(z_space,r_space);

% psi = get_psi(B_z,r_probe,t+4);
% psi = get_psi(B_z,r_probe,t+2);
psi = get_psi(B_z,r_probe,t);
psi = griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r,'cubic');

psi_min = min(min(psi));
psi_max = max(max(psi));
contour_layer = linspace(psi_min,psi_max,20);

% if show_xpoint
%     [~,~,pos_xz,pos_xr,~,~] = search_xo(psi,z_space,r_space);
%     dz = 0.02;
%     dr = 0.03;
%     pos_xz_lower = pos_xz - dz;
%     pos_xr_lower = pos_xr - dr;
% end

f = figure;
f.Units = 'normalized';
f.Position = [0.1,0.2,0.8,0.4];
pos1 = [0.07,0.2,0.35,0.6];
pos2 = [0.58,0.2,0.35,0.6];

subplot('Position',pos1);
% subplot(1,2,1);
[SXR_mesh_z1,SXR_mesh_r1] = meshgrid(z_space_SXR1,r_space_SXR1);
[~,h1] = contourf(SXR_mesh_z1,SXR_mesh_r1,EE1,20);
h1.LineStyle = 'none';
clim([0,0.2]);
% clim([0,0.4]);
% c=colorbar('westoutside');
c=colorbar;
c.Label.String='Intensity [a.u.]';c.FontSize=18;
hold on

[~,hp1]=contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'-k','Fill','off');
% if layer
    hp1.LineWidth = 1.5;
% end

if show_localmax
    localmax_idx = imregionalmax(EE1);
%     localmax_pos_r = SXR_mesh_r1(localmax_idx);
%     localmax_pos_z = SXR_mesh_z1(localmax_idx);
    EE_localmax = EE1.*localmax_idx;
    [~, localmax_idx] = maxk(EE_localmax(:),2);
    localmax_pos_r = SXR_mesh_r1(localmax_idx);
    localmax_pos_z = SXR_mesh_z1(localmax_idx);
    plot(localmax_pos_z,localmax_pos_r,'r*');
end

if show_xpoint
    [~,~,pos_xz,pos_xr,~,~] = search_xo(psi,z_space,r_space);
    dz = 0.02;
    dr = 0.03;
    pos_xz_lower = pos_xz - dz;
    pos_xr_lower = pos_xr - dr;
    r = rectangle('Position',[pos_xz_lower pos_xr_lower dz*2 dr*2]);
    r.EdgeColor = 'red';
    r.LineWidth = 1.5;
end
title(strcat(num2str(t),' us'));
xlabel('z [m]');
ylabel('r [m]');
ax = gca;
ax.FontSize = 18; 
hold off

subplot('Position',pos2);
% subplot(1,2,2);
[SXR_mesh_z2,SXR_mesh_r2] = meshgrid(z_space_SXR2,r_space_SXR2);
[~,h2] = contourf(SXR_mesh_z2,SXR_mesh_r2,EE2,20);
h2.LineStyle = 'none';
% clim([0,0.15]);
% clim([0,0.2]);
clim([0,0.1]);
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
hold on

[~,hp2]=contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'-k','Fill','off');
% if layer
    hp2.LineWidth = 1.5;
% end

if show_xpoint
    r = rectangle('Position',[pos_xz_lower pos_xr_lower dz*2 dr*2]);
    r.EdgeColor = 'red';
    r.LineWidth = 1.5;
end
title(strcat(num2str(t),' us'));
xlabel('z [m]');
ylabel('r [m]');
ax = gca;
ax.FontSize = 18; 
hold off

if save
    pathname = '/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/ReconstructionResults/';
    if filter & NL
        directory = '/NLF_NLR/';
    elseif ~filter & NL
        directory = '/LF_NLR/';
    elseif filter & ~NL
        directory = '/NLF_LR/';
    else
        directory = '/LF_LR/';
    end
    foldername = strcat(pathname,directory,num2str(date),'/shot',num2str(shot),'_wide');
    if exist(foldername,'dir') == 0
        mkdir(foldername);
    end
    filename = strcat('/shot',num2str(shot),'_',num2str(t),'us.png');
    saveas(gcf,strcat(foldername,filename));
end
end