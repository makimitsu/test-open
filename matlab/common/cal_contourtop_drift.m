function [Contopdata] = cal_contourtop_drift(CONTOP,PCB,FIG,pathname,savename,range)

savename_contop = [pathname.mat,'/contourtop/',num2str(PCB.date),'_a039_',num2str(PCB.shot),'_',num2str(FIG.start),'_',num2str(FIG.dt),'_',num2str(FIG.tate*FIG.yoko),'.mat'];
if exist(savename_contop,"file")
    load(savename_contop,'Contopdata')
else
    if exist(savename.highmesh_psi,"file")
        load(savename.highmesh_psi,'highmesh_PCBgrid2D','highmesh_PCBdata2D')
    else
        warning([savename.highmesh_psi, 'does not exist.'])
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
    low_z = downsample(z,CONTOP.ds_rate);
    low_r = downsample(r,CONTOP.ds_rate);
    Contopdata.rq = low_r;
    Contopdata.zq = zeros(size(Contopdata.rq,1),FIG.tate*FIG.yoko);
    Contopdata.trange = zeros(FIG.tate*FIG.yoko,1);
    Contopdata.Vcontop_r = zeros(size(Contopdata.rq,1),FIG.tate*FIG.yoko);%磁力線移動速度のr成分[km/s]

    for i_t = 1:FIG.tate*FIG.yoko
        pcb_time = FIG.start+FIG.dt*(i_t-1);
        Contopdata.trange(i_t,1) = pcb_time;
        idx_highmesh_t = knnsearch(highmesh_PCBdata2D.trange',pcb_time);
        psi = highmesh_PCBdata2D.psi(:,:,idx_highmesh_t);
        psi = psi(idx_r_min:idx_r_cal_max,idx_z_min:idx_z_max);
        low_psi = downsample(psi',CONTOP.ds_rate);
        low_psi = downsample(low_psi',CONTOP.ds_rate);
        %X点を探索
        for j = 1:size(low_psi,2)
            if j == 1
                [saddle_psi,saddle_r] = max(low_psi(:,j));
                saddle_z = j;
            else
                if max(low_psi(:,j)) < saddle_psi
                    [saddle_psi,saddle_r] = max(low_psi(:,j));
                    saddle_z = j;
                end
            end
        end
        psi_before = highmesh_PCBdata2D.psi(:,:,idx_highmesh_t-1);%1us前のpsi
        psi_before = psi_before(idx_r_min:idx_r_cal_max,idx_z_min:idx_z_max);
        psi_after = highmesh_PCBdata2D.psi(:,:,idx_highmesh_t+1);%1us後のpsi
        psi_after = psi_after(idx_r_min:idx_r_cal_max,idx_z_min:idx_z_max);
        for i_r = 1:size(low_r,1)
            [level,idx_z] = min(low_psi(i_r,:));%追跡する磁力線のpsi
            Contopdata.zq(i_r,i_t) = low_z(idx_z);
            %1us前の同磁力線位置を探索
            M_before = contourc(z,r,psi_before,[level,level]);
            IDX_start = find(M_before(1,:)==level);
            n_line = size(IDX_start,2);%同じlevelの連続曲線数
            contop_r_before = zeros(n_line,1);
            for i = 1:n_line
                if i < n_line
                    M_line = M_before(:,IDX_start(i)+1:IDX_start(i+1)-1);
                else
                    M_line = M_before(:,IDX_start(i)+1:end);
                end
                ext_r_before = [max(M_line(2,:));min(M_line(2,:))];
                idx_contop_r_before = knnsearch(ext_r_before,low_r(saddle_r));
                contop_r_before(i,1) = ext_r_before(idx_contop_r_before);
            end
            idx_r_before = knnsearch(contop_r_before,low_r(i_r));
            r_before = contop_r_before(idx_r_before,1);
            %1us後の同磁力線位置を探索
            M_after = contourc(z,r,psi_after,[level,level]);
            IDX_start = find(M_after(1,:)==level);
            n_line = size(IDX_start,2);%同じlevelの連続曲線数
            contop_r_after = zeros(n_line,1);
            for i = 1:n_line
                if i < n_line
                    M_line = M_after(:,IDX_start(i)+1:IDX_start(i+1)-1);
                else
                    M_line = M_after(:,IDX_start(i)+1:end);
                end
                ext_r_after = [max(M_line(2,:));min(M_line(2,:))];
                idx_contop_r_after = knnsearch(ext_r_after,low_r(saddle_r));
                contop_r_after(i,1) = ext_r_after(idx_contop_r_after);
            end
            idx_r_after = knnsearch(contop_r_after,low_r(i_r));
            r_after = contop_r_after(idx_r_after,1);
            if (size(r_before,1) == 1) && (size(r_after,1) == 1)
                Contopdata.Vcontop_r(i_r,i_t) = (r_after-r_before)/(2E-6)*1E-3;
            else
                warning('stop')
                return
            end
        end
        Contopdata.Vcontop_r(:,i_t) = filloutliers(Contopdata.Vcontop_r(:,i_t),"linear");
    end
    idx_r_max = knnsearch(Contopdata.rq(:,1),range.r_max);
    Contopdata.rq(idx_r_max:end,:) = [];
    Contopdata.zq(idx_r_max:end,:) = [];
    Contopdata.Vcontop_r(idx_r_max:end,:,:) = [];
    save(savename_contop,'Contopdata')
end
