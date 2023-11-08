function [Curvedata2D] = cal_curve_drift(CURVE,PCB,FIG,pathname,savename,range)

savename_curve = [pathname.mat,'/curve/',num2str(PCB.date),'_a039_',num2str(PCB.shot),'_',num2str(FIG.start),'_',num2str(FIG.dt),'_',num2str(FIG.tate*FIG.yoko),'.mat'];
if exist(savename_curve,"file")
    load(savename_curve,'Curvedata2D')
else
    if exist(savename.highmesh_psi,"file")
        load(savename.highmesh_psi,'highmesh_PCBgrid2D','highmesh_PCBdata2D')
    else
        warning([savename.highmesh_psi, 'does not exist.'])
        return
    end
    if exist(savename.pcb,"file")
        load(savename.pcb,'PCBgrid2D','PCBdata2D')
    else
        warning([savename.pcb, 'does not exist.'])
        return
    end
    m_i = 1.67E-27*CURVE.A;%イオン質量[kg]
    q_i = 1.6E-19;%イオン電荷[C]
    z = highmesh_PCBgrid2D.zq(1,:);
    r = highmesh_PCBgrid2D.rq(:,1);
    idx_z_max = knnsearch(z',range.z_max);
    idx_z_min = knnsearch(z',range.z_min);
    idx_r_max = knnsearch(r,range.r_max);
    idx_r_min = knnsearch(r,range.r_min);
    z = z(idx_z_min:idx_z_max)';
    r = r(idx_r_min:idx_r_max);
    low_z = downsample(z,CURVE.ds_rate);
    low_r = downsample(r,CURVE.ds_rate);
    [Curvedata2D.zq,Curvedata2D.rq] = meshgrid(low_z,low_r);
    Curvedata2D.trange = zeros(FIG.tate*FIG.yoko,1);
    Curvedata2D.Fcurve_z = zeros(size(Curvedata2D.zq,1),size(Curvedata2D.zq,2),FIG.tate*FIG.yoko);
    Curvedata2D.Fcurve_r = zeros(size(Curvedata2D.zq,1),size(Curvedata2D.zq,2),FIG.tate*FIG.yoko);
    Curvedata2D.Vcurve_z = zeros(size(Curvedata2D.zq,1),size(Curvedata2D.zq,2),FIG.tate*FIG.yoko);%曲率ドリフトのz成分[km/s]
    Curvedata2D.Vcurve_r = zeros(size(Curvedata2D.zq,1),size(Curvedata2D.zq,2),FIG.tate*FIG.yoko);%曲率ドリフトのr成分[km/s]
    Curvedata2D.Vgc = zeros(size(Curvedata2D.zq,1),size(Curvedata2D.zq,2),FIG.tate*FIG.yoko);%[m/s]
    Curvedata2D.Vgc_p = zeros(size(Curvedata2D.zq,1),size(Curvedata2D.zq,2),FIG.tate*FIG.yoko);%[m/s]
    Curvedata2D.Vgc_t = zeros(size(Curvedata2D.zq,1),size(Curvedata2D.zq,2),FIG.tate*FIG.yoko);%[m/s]

    for i_t = 1:FIG.tate*FIG.yoko
        pcb_time = FIG.start+FIG.dt*(i_t-1);
        Curvedata2D.trange(i_t,1) = pcb_time;
        idx_pcb_t = knnsearch(PCBdata2D.trange',pcb_time);
        idx_highmesh_t = knnsearch(highmesh_PCBdata2D.trange',pcb_time);
        psi = highmesh_PCBdata2D.psi(:,:,idx_highmesh_t);
        psi = psi(idx_r_min:idx_r_max,idx_z_min:idx_z_max);
        low_psi = downsample(psi',CURVE.ds_rate);
        low_psi = downsample(low_psi',CURVE.ds_rate);
        for i_r = 1:size(low_r,1)
            for i_z = 1:size(low_z,1)
                idx_pcb_r = knnsearch(PCBgrid2D.rq(:,1),low_r(i_r));
                idx_pcb_z = knnsearch(PCBgrid2D.zq(1,:)',low_z(i_z));
                level = low_psi(i_r,i_z);
                M = contourc(z,r,psi,[level,level]);
                IDX_start = find(M(1,:)==level);
                n_line = size(IDX_start,2);%同じlevelの連続曲線数
                for i = 1:n_line
                    if i < n_line
                        M_line = M(:,IDX_start(i)+1:IDX_start(i+1)-1);
                    else
                        M_line = M(:,IDX_start(i)+1:end);
                    end
                    IDX_center = find((M_line(1,:)==low_z(i_z)) & (M_line(2,:)==low_r(i_r)));
                    if isempty(IDX_center)
                        IDX_center = 0;
                    elseif size(IDX_center,2)>1
                        IDX_center = IDX_center(1,1);
                    end
                    if (IDX_center-CURVE.c_width > 1) && (IDX_center+CURVE.c_width < size(M_line,2))
                        M_center = M_line(:,IDX_center-CURVE.c_width:IDX_center+CURVE.c_width);
                        if size(M_center,2) == 2*CURVE.c_width+1
                            z_c = M_center(1,:);
                            r_c = M_center(2,:);
                            X = [z_c',r_c'];
                            [~,~,K2] = curvature(X);
                            sum_Etdt = sum(PCBdata2D.Et(idx_pcb_r,idx_pcb_z,idx_pcb_t-10:idx_pcb_t),3)*1E-6;
                            Curvedata2D.Vgc_t(i_r,i_z,i_t) = q_i/m_i*sum_Etdt;%Toroidal component of Gyrocenter velocity[m/s]
                            Bt_ext = PCBdata2D.Bt_ext(idx_pcb_r,idx_pcb_z,idx_pcb_t);
                            Bt = PCBdata2D.Bt(idx_pcb_r,idx_pcb_z,idx_pcb_t);
                            Br = PCBdata2D.Br(idx_pcb_r,idx_pcb_z,idx_pcb_t);
                            Bz = PCBdata2D.Bz(idx_pcb_r,idx_pcb_z,idx_pcb_t);
                            Bp = sqrt(Br^2+Bz^2);
                            % absB = sqrt(Bt_ext^2+Bp^2);
                            absB = sqrt(Bt^2+Bp^2);
                            % Curvedata2D.Vgc_p(i_r,i_z,i_t) = Curvedata2D.Vgc_t(i_r,i_z,i_t)*Bp/Bt_ext;%Poloidal component of Gyrocenter velocity[m/s]
                            Curvedata2D.Vgc_p(i_r,i_z,i_t) = Curvedata2D.Vgc_t(i_r,i_z,i_t)*Bp/Bt;%Poloidal component of Gyrocenter velocity[m/s]
                            Fcurve_z = -K2(:,1)*m_i*Curvedata2D.Vgc_p(i_r,i_z,i_t)^2;
                            Fcurve_z = filloutliers(Fcurve_z,"linear");
                            Fcurve_r = -K2(:,2)*m_i*Curvedata2D.Vgc_p(i_r,i_z,i_t)^2;
                            Fcurve_r = filloutliers(Fcurve_r,"linear");
                            Curvedata2D.Fcurve_z(i_r,i_z,i_t) = Fcurve_z(CURVE.c_width+1);%遠心力のz成分[N]
                            Curvedata2D.Fcurve_r(i_r,i_z,i_t) = Fcurve_r(CURVE.c_width+1)+m_i/low_r(i_r)*Curvedata2D.Vgc_t(i_r,i_z,i_t)^2;%遠心力のr成分[N]
                            % Curvedata2D.Vcurve_z(i_r,i_z,i_t) = Curvedata2D.Fcurve_r(i_r,i_z,i_t)*Bt_ext/(q_i*absB^2)*1E-3;
                            % Curvedata2D.Vcurve_r(i_r,i_z,i_t) = -Curvedata2D.Fcurve_z(i_r,i_z,i_t)*Bt_ext/(q_i*absB^2)*1E-3;
                            Curvedata2D.Vcurve_z(i_r,i_z,i_t) = Curvedata2D.Fcurve_r(i_r,i_z,i_t)*Bt/(q_i*absB^2)*1E-3;
                            Curvedata2D.Vcurve_r(i_r,i_z,i_t) = -Curvedata2D.Fcurve_z(i_r,i_z,i_t)*Bt/(q_i*absB^2)*1E-3;
                        end
                    end
                end
            end
        end
    end
    Curvedata2D.zq([1 end],:,:) = [];
    Curvedata2D.zq(:,[1 end],:) = [];
    Curvedata2D.rq([1 end],:,:) = [];
    Curvedata2D.rq(:,[1 end],:) = [];
    Curvedata2D.Fcurve_z([1 end],:,:) = [];
    Curvedata2D.Fcurve_z(:,[1 end],:) = [];
    Curvedata2D.Fcurve_r([1 end],:,:) = [];
    Curvedata2D.Fcurve_r(:,[1 end],:) = [];
    Curvedata2D.Vcurve_z([1 end],:,:) = [];
    Curvedata2D.Vcurve_z(:,[1 end],:) = [];
    Curvedata2D.Vcurve_r([1 end],:,:) = [];
    Curvedata2D.Vcurve_r(:,[1 end],:) = [];
    Curvedata2D.Vgc([1 end],:,:) = [];
    Curvedata2D.Vgc(:,[1 end],:) = [];
    Curvedata2D.Vgc_p([1 end],:,:) = [];
    Curvedata2D.Vgc_p(:,[1 end],:) = [];
    Curvedata2D.Vgc_t([1 end],:,:) = [];
    Curvedata2D.Vgc_t(:,[1 end],:) = [];
    Curvedata2D.Fcurve = sqrt(Curvedata2D.Fcurve_z.^2 + Curvedata2D.Fcurve_r.^2);
    Curvedata2D.Vcurve = sqrt(Curvedata2D.Vcurve_z.^2 + Curvedata2D.Vcurve_r.^2);
    save(savename_curve,'Curvedata2D')
end