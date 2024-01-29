function [Maglinedata2D] = cal_magline_drift(MAGLINE,PCB,FIG,pathname,savename,range)

savename_magline = [pathname.mat,'/magline/',num2str(PCB.date),'_a039_',num2str(PCB.shot),'_',num2str(FIG.start),'_',num2str(FIG.dt),'_',num2str(FIG.tate*FIG.yoko),'.mat'];
if exist(savename_magline,"file")
    load(savename_magline,'Maglinedata2D')
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
    z = highmesh_PCBgrid2D.zq(1,:);
    r = highmesh_PCBgrid2D.rq(:,1);
    idx_z_max = knnsearch(z',range.z_max);
    idx_z_min = knnsearch(z',range.z_min);
    idx_r_cal_max = knnsearch(r,range.r_cal_max);
    idx_r_min = knnsearch(r,range.r_min);
    z = z(idx_z_min:idx_z_max)';
    r = r(idx_r_min:idx_r_cal_max);
    low_z = downsample(z,MAGLINE.ds_rate);
    low_r = downsample(r,MAGLINE.ds_rate);
    [Maglinedata2D.zq,Maglinedata2D.rq] = meshgrid(low_z,low_r);
    Maglinedata2D.trange = zeros(FIG.tate*FIG.yoko,1);
    Maglinedata2D.Vmagline_z = zeros(size(Maglinedata2D.zq,1),size(Maglinedata2D.zq,2),FIG.tate*FIG.yoko);%磁力線移動速度のz成分[km/s] (=0と仮定)
    Maglinedata2D.Vmagline_r = zeros(size(Maglinedata2D.zq,1),size(Maglinedata2D.zq,2),FIG.tate*FIG.yoko);%磁力線移動速度のr成分[km/s]

    for i_t = 1:FIG.tate*FIG.yoko
        pcb_time = FIG.start+FIG.dt*(i_t-1);
        Maglinedata2D.trange(i_t,1) = pcb_time;
        idx_highmesh_t = knnsearch(highmesh_PCBdata2D.trange',pcb_time);
        psi = highmesh_PCBdata2D.psi(:,:,idx_highmesh_t);
        psi = psi(idx_r_min:idx_r_cal_max,idx_z_min:idx_z_max);
        low_psi = downsample(psi',MAGLINE.ds_rate);
        low_psi = downsample(low_psi',MAGLINE.ds_rate);
        psi_before = highmesh_PCBdata2D.psi(:,:,idx_highmesh_t-1);%1us前のpsi
        psi_before = psi_before(idx_r_min:idx_r_cal_max,idx_z_min:idx_z_max);
        psi_after = highmesh_PCBdata2D.psi(:,:,idx_highmesh_t+1);%1us後のpsi
        psi_after = psi_after(idx_r_min:idx_r_cal_max,idx_z_min:idx_z_max);
        for i_r = 1:size(low_r,1)
            for i_z = 1:size(low_z,1)
                level = low_psi(i_r,i_z);%追跡する磁力線のpsi
                %1us前の同磁力線位置を探索
                M_before = contourc(z,r,psi_before,[level,level]);
                IDX_start = find(M_before(1,:)==level);
                n_line = size(IDX_start,2);%同じlevelの連続曲線数
                magline_r_before = zeros(n_line,1);
                for i = 1:n_line
                    if i < n_line
                        M_line = M_before(:,IDX_start(i)+1:IDX_start(i+1)-1);
                    else
                        M_line = M_before(:,IDX_start(i)+1:end);
                    end
                    idx_before = knnsearch(M_line(1,:)',low_z(i_z));
                    magline_r_before(i,1) = M_line(2,idx_before);
                end
                idx_r_before = knnsearch(magline_r_before,low_r(i_r));
                r_before = magline_r_before(idx_r_before,1);
                %1us後の同磁力線位置を探索
                M_after = contourc(z,r,psi_after,[level,level]);
                IDX_start = find(M_after(1,:)==level);
                n_line = size(IDX_start,2);%同じlevelの連続曲線数
                magline_r_after = zeros(n_line,1);
                for i = 1:n_line
                    if i < n_line
                        M_line = M_after(:,IDX_start(i)+1:IDX_start(i+1)-1);
                    else
                        M_line = M_after(:,IDX_start(i)+1:end);
                    end
                    idx_after = knnsearch(M_line(1,:)',low_z(i_z));
                    magline_r_after(i,1) = M_line(2,idx_after);
                end
                idx_r_after = knnsearch(magline_r_after,low_r(i_r));
                r_after = magline_r_after(idx_r_after,1);
                if (size(r_before,1) == 1) && (size(r_after,1) == 1)
                    Maglinedata2D.Vmagline_r(i_r,i_z,i_t) = (r_after-r_before)/(2E-6)*1E-3;
                end
            end
        end
    end
    idx_r_max = knnsearch(Maglinedata2D.rq(:,1),range.r_max);
    Maglinedata2D.zq(idx_r_max:end,:) = [];
    Maglinedata2D.rq(idx_r_max:end,:) = [];
    Maglinedata2D.Vmagline_z(idx_r_max:end,:,:) = [];
    Maglinedata2D.Vmagline_r(idx_r_max:end,:,:) = [];
    save(savename_magline,'Maglinedata2D')
end