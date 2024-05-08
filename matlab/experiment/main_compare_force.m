clear all
addpath '/Users/rsomeya/Documents/lab/matlab/common';
addpath '/Users/rsomeya/Documents/lab/matlab/common'/curvature;
run define_path.m

PCB.date = 230830;

CURVE.c_width = 5;
CURVE.ds_rate = 5;
CURVE.A = 40;
NABLAB.A = 40;
NABLAB.T_i = 10;
% range.z_min = 0;
% range.z_max = 0.05;
range.z_min = -0.05;
range.z_max = 0.1;
range.r_min = 0.11;
range.r_max = 0.27;
FIG.tate = 1;%【input】プロット枚数(縦)
FIG.yoko = 1;%【input】プロット枚数(横)
FIG.dt = 1;%【input】プロット時間間隔[us]
cal_type = 'curve';%'curve','ExB','nablaB'
plot_type = 'F';%'V_zr','V_z','V_r','V','F_zr','F_z','F_r','F'

switch PCB.date
    case 230828
        PCB.shot = 2301;
        FIG.start = 482;
        savename.highmesh_psi = [pathname.mat,'/pcb_processed/','mesh500_a039_2301.mat'];
        savename.pcb = [pathname.mat,'/pcb_processed/','a039_2301.mat'];
        savename.ESP = [pathname.mat,'/ESP/','mesh21_230828_shot5-63.mat'];
    case 230830
        PCB.shot = 2437;
        FIG.start = 470;
        savename.highmesh_psi = [pathname.mat,'/pcb_processed/','mesh500_a039_2437.mat'];
        savename.pcb = [pathname.mat,'/pcb_processed/','a039_2437.mat'];
        savename.ESP = [pathname.mat,'/ESP/','mesh21_230830_shot11-60.mat'];
    case 230914
        PCB.shot = 2862;
        FIG.start = 500;
        savename.highmesh_psi = [pathname.mat,'/pcb_processed/','mesh500_a039_2862.mat'];
        savename.pcb = [pathname.mat,'/pcb_processed/','a039_2862.mat'];
    otherwise
        warning('PCB.date Error')
        return
end


switch cal_type
    case 'curve'
        [Curvedata2D] = cal_curve_drift(CURVE,PCB,FIG,pathname,savename,range);
        plot_curve_drift(Curvedata2D,plot_type,FIG,savename)
    case 'nablaB'
        [NablaBdata2D] = cal_nablaB_drift(NABLAB,PCB,FIG,pathname,savename,range);
        plot_nablaB_drift(NablaBdata2D,plot_type,FIG,savename)
    case 'ExB'
        if exist(savename.ESP,"file")
            load(savename.ESP,'ESPdata2D')
        else
            warning([savename.ESP, 'does not exist.'])
            return
        end
        idx_time = knnsearch(ESPdata2D.trange',FIG.start);
        z = ESPdata2D.zq(1,:);
        r = ESPdata2D.rq(:,1);
        idx_range.z_max = knnsearch(z',range.z_max);
        idx_range.z_min = knnsearch(z',range.z_min);
        idx_range.r_max = knnsearch(r,range.r_max);
        idx_range.r_min = knnsearch(r,range.r_min);
        z1 = z(idx_range.z_min:idx_range.z_max)';
        r1 = r(idx_range.r_min:idx_range.r_max);
        [mesh_z,mesh_r] = meshgrid(z1,r1);
        E_z = squeeze(ESPdata2D.Ez(idx_time,idx_range.r_min:idx_range.r_max,idx_range.z_min:idx_range.z_max));
        E_r = squeeze(ESPdata2D.Er(idx_time,idx_range.r_min:idx_range.r_max,idx_range.z_min:idx_range.z_max));
        absE = sqrt(E_z.^2+E_r.^2);
        Fcurve_z = q_i*E_z;
        Fcurve_r = q_i*E_r;
        Fcurve = q_i*absE;
        figure('Position',[0 0 1000*(range.z_max-range.z_min)/(range.r_max-range.r_min)+100 1500],'visible','on')
        switch plot_type
            case 'F_r'
                contourf(mesh_z,mesh_r,Fcurve_r,100,'edgecolor','none')
                c = colorbar;
                colormap(jet)
                c.Label.String = 'R component of Electrical Force [N]';
            case 'F_z'
                contourf(mesh_z,mesh_r,Fcurve_z,100,'edgecolor','none')
                c = colorbar;
                colormap(jet)
                clim([-2E-16 2E-16])
                c.Label.String = 'Z component of Electrical Force [N]';
            case 'F'
                contourf(mesh_z,mesh_r,Fcurve,100,'edgecolor','none')
                c = colorbar;
                colormap(jet)
                clim([0 3E-16])
                c.Label.String = 'Strength of Electrical Force [N]';
        end
end

hold on

% if exist(savename.highmesh_psi,"file")
%     load(savename.highmesh_psi,'highmesh_PCBgrid2D','highmesh_PCBdata2D')
%     idx_pcb = knnsearch(highmesh_PCBdata2D.trange',FIG.start);
%     contour(highmesh_PCBgrid2D.zq(1,:),highmesh_PCBgrid2D.rq(:,1),squeeze(highmesh_PCBdata2D.psi(:,:,idx_pcb)),[-20e-3:0.05e-3:40e-3],'black','LineWidth',1)
% else
%     warning([savename.highmesh_psi, 'does not exist.'])
%     return
% end

% daspect([1 1 1])
% xlim([range.z_min range.z_max])
% ylim([range.r_min range.r_max])
% xlabel('Z [m]')
% ylabel('R [m]')
% ax = gca;
% ax.FontSize = 18;
