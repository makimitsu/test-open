function plot_save_sxr(grid2D,data2D,range,date,shot,t,EE1,EE2,EE3,EE4,show_localmax,show_xpoint,save,filter,NL)

f = figure;
f.Units = 'normalized';
f.Position = [0.1,0.2,0.8,0.8];

range = range./1000;
zmin = range(1);
zmax = range(2);
rmin = range(5);
rmax = range(6);
r_space_SXR = linspace(rmin,rmax,size(EE1,1));
z_space_SXR = linspace(zmin,zmax,size(EE1,2));

r_range = find(0.060<=r_space_SXR & r_space_SXR<=0.330);
r_space_SXR = r_space_SXR(r_range);
z_range = find(-0.17<=z_space_SXR & z_space_SXR<=0.17);

z_space_SXR = z_space_SXR(z_range);

EE1 = EE1(r_range,z_range);
EE2 = EE2(r_range,z_range);
EE3 = EE3(r_range,z_range);
EE4 = EE4(r_range,z_range);

psi_mesh_z = grid2D.zq;
psi_mesh_r = grid2D.rq;
t_idx = find(data2D.trange==t);
psi = data2D.psi(:,:,t_idx);

psi_min = min(min(psi));
psi_max = max(max(psi));
contour_layer = linspace(psi_min,psi_max,20);

subplot(2,2,1);
[SXR_mesh_z,SXR_mesh_r] = meshgrid(z_space_SXR,r_space_SXR);
[~,h1] = contourf(SXR_mesh_z,SXR_mesh_r,EE1,20);
h1.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
hold on
[~,hp1]=contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'-k','Fill','off');
hp1.LineWidth = 1.5;
if show_localmax
    localmax_idx = imregionalmax(EE1);
    EE_localmax = EE1.*localmax_idx;
    [~, localmax_idx] = maxk(EE_localmax(:),2);
    localmax_pos_r = SXR_mesh_r(localmax_idx);
    localmax_pos_z = SXR_mesh_z(localmax_idx);
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
title('1');

subplot(2,2,2);
[~,h2] = contourf(SXR_mesh_z,SXR_mesh_r,EE2,20);
h2.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
hold on
[~,hp2]=contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'-k','Fill','off');
hp2.LineWidth = 1.5;
if show_localmax
    localmax_idx = imregionalmax(EE2);
    EE_localmax = EE2.*localmax_idx;
    [~, localmax_idx] = maxk(EE_localmax(:),2);
    localmax_pos_r = SXR_mesh_r(localmax_idx);
    localmax_pos_z = SXR_mesh_z(localmax_idx);
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
title('2');

subplot(2,2,3);
[~,h3] = contourf(SXR_mesh_z,SXR_mesh_r,EE3,20);
h3.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
hold on
[~,hp3]=contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'-k','Fill','off');
hp3.LineWidth = 1.5;
if show_localmax
    localmax_idx = imregionalmax(EE3);
    EE_localmax = EE3.*localmax_idx;
    [~, localmax_idx] = maxk(EE_localmax(:),2);
    localmax_pos_r = SXR_mesh_r(localmax_idx);
    localmax_pos_z = SXR_mesh_z(localmax_idx);
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
title('3');

subplot(2,2,4);
[~,h4] = contourf(SXR_mesh_z,SXR_mesh_r,EE4,20);
h4.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
hold on
[~,hp4]=contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'-k','Fill','off');
hp4.LineWidth = 1.5;
if show_localmax
    localmax_idx = imregionalmax(EE4);
    EE_localmax = EE4.*localmax_idx;
    [~, localmax_idx] = maxk(EE_localmax(:),2);
    localmax_pos_r = SXR_mesh_r(localmax_idx);
    localmax_pos_z = SXR_mesh_z(localmax_idx);
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
title('4');

if save
    pathname = getenv('SXR_RECONSTRUCTED_DIR');
    if filter & NL
        directory = '/NLF_NLR/';
    elseif ~filter & NL
        directory = '/LF_NLR/';
    elseif filter & ~NL
        directory = '/NLF_LR/';
    else
        directory = '/LF_LR/';
    end
    foldername = strcat(pathname,directory,'/',num2str(date),'/shot',num2str(shot),'_wide');
    if exist(foldername,'dir') == 0
        mkdir(foldername);
    end
    filename = strcat('/shot',num2str(shot),'_',num2str(t),'us.png');
    saveas(gcf,strcat(foldername,filename));
end

end