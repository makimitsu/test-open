function [NablaBdata2D] = cal_nablaB_drift(NABLAB,PCB,FIG,pathname,savename,range)

savename_nablaB = [pathname.mat,'/nablaB/',num2str(PCB.date),'_a039_',num2str(PCB.shot),'_',num2str(FIG.start),'_',num2str(FIG.dt),'_',num2str(FIG.tate*FIG.yoko),'.mat'];
if exist(savename_nablaB,"file")
    load(savename_nablaB,'NablaBdata2D')
else
    if exist(savename.pcb,"file")
        load(savename.pcb,'PCBgrid2D','PCBdata2D')
    else
        warning([savename.pcb, 'does not exist.'])
        return
    end
    m_i = 1.67E-27*NABLAB.A;%イオン質量[kg]
    q_i = 1.6E-19;%イオン電荷[C]
    V_L = 1.4E4*sqrt(NABLAB.T_i/NABLAB.A);%Gyration velocity[m/s]
    NablaBdata2D.V_L = V_L;
    z = PCBgrid2D.zq(1,:);
    r = PCBgrid2D.rq(:,1);
    idx_z_max = knnsearch(z',range.z_max);
    idx_z_min = knnsearch(z',range.z_min);
    idx_r_max = knnsearch(r,range.r_max);
    idx_r_min = knnsearch(r,range.r_min);
    z = z(idx_z_min:idx_z_max)';
    r = r(idx_r_min:idx_r_max);
    [NablaBdata2D.zq,NablaBdata2D.rq] = meshgrid(z,r);
    NablaBdata2D.trange = zeros(FIG.tate*FIG.yoko,1);
    NablaBdata2D.FnablaB_z = zeros(size(NablaBdata2D.zq,1),size(NablaBdata2D.zq,2),FIG.tate*FIG.yoko);
    NablaBdata2D.FnablaB_r = zeros(size(NablaBdata2D.zq,1),size(NablaBdata2D.zq,2),FIG.tate*FIG.yoko);
    NablaBdata2D.FnablaB = zeros(size(NablaBdata2D.zq,1),size(NablaBdata2D.zq,2),FIG.tate*FIG.yoko);
    NablaBdata2D.VnablaB_z = zeros(size(NablaBdata2D.zq,1),size(NablaBdata2D.zq,2),FIG.tate*FIG.yoko);
    NablaBdata2D.VnablaB_r = zeros(size(NablaBdata2D.zq,1),size(NablaBdata2D.zq,2),FIG.tate*FIG.yoko);
    NablaBdata2D.VnablaB = zeros(size(NablaBdata2D.zq,1),size(NablaBdata2D.zq,2),FIG.tate*FIG.yoko);
    dz = PCBgrid2D.zq(1,2) - PCBgrid2D.zq(1,1);
    dr = PCBgrid2D.rq(2,1) - PCBgrid2D.rq(1,1);
    for i_t = 1:FIG.tate*FIG.yoko
        pcb_time = FIG.start+FIG.dt*(i_t-1);
        NablaBdata2D.trange(i_t,1) = pcb_time;
        idx_pcb_t = knnsearch(PCBdata2D.trange',pcb_time);
        Bz = PCBdata2D.Bz(:,:,idx_pcb_t);
        Br = PCBdata2D.Br(:,:,idx_pcb_t);
        Bt_ext = PCBdata2D.Bt_ext(:,:,idx_pcb_t);
        Bt = PCBdata2D.Bt(:,:,idx_pcb_t);
        Bz = Bz(idx_r_min:idx_r_max,idx_z_min:idx_z_max);
        Br = Br(idx_r_min:idx_r_max,idx_z_min:idx_z_max);
        Bt_ext = Bt_ext(idx_r_min:idx_r_max,idx_z_min:idx_z_max);
        Bt = Bt(idx_r_min:idx_r_max,idx_z_min:idx_z_max);
        % absB = sqrt(Bz.^2+Br.^2+Bt_ext.^2);
        absB = sqrt(Bz.^2+Br.^2+Bt.^2);
        NablaBdata2D.FnablaB_z(:,1:end-1,i_t) = -m_i*V_L^2./(2*absB(:,1:end-1)).*diff(absB(:,:),1,2)./dz;%∇B力のz成分[N]
        NablaBdata2D.FnablaB_r(1:end-1,:,i_t) = -m_i*V_L^2./(2*absB(1:end-1,:)).*diff(absB(:,:),1,1)./dr;%∇B力のr成分[N]
        NablaBdata2D.FnablaB(:,:,i_t) = sqrt(NablaBdata2D.FnablaB_z(:,:,i_t).^2 + NablaBdata2D.FnablaB_r(:,:,i_t).^2);
        % NablaBdata2D.VnablaB_z(:,:,i_t) = NablaBdata2D.FnablaB_r(:,:,i_t).*Bt_ext./(q_i*absB.^2)*1E-3;%∇Bドリフトのz成分[km/s]
        % NablaBdata2D.VnablaB_r(:,:,i_t) = -NablaBdata2D.FnablaB_z(:,:,i_t).*Bt_ext./(q_i*absB.^2)*1E-3;%∇Bドリフトのr成分[km/s]
        NablaBdata2D.VnablaB_z(:,:,i_t) = NablaBdata2D.FnablaB_r(:,:,i_t).*Bt./(q_i*absB.^2)*1E-3;%∇Bドリフトのz成分[km/s]
        NablaBdata2D.VnablaB_r(:,:,i_t) = -NablaBdata2D.FnablaB_z(:,:,i_t).*Bt./(q_i*absB.^2)*1E-3;%∇Bドリフトのr成分[km/s]
        NablaBdata2D.VnablaB(:,:,i_t) = sqrt(NablaBdata2D.VnablaB_z(:,:,i_t).^2 + NablaBdata2D.VnablaB_r(:,:,i_t).^2);
    end
    NablaBdata2D.zq(end,:,:) = [];
    NablaBdata2D.zq(:,end,:) = [];
    NablaBdata2D.rq(end,:,:) = [];
    NablaBdata2D.rq(:,end,:) = [];
    NablaBdata2D.FnablaB_z(end,:,:) = [];
    NablaBdata2D.FnablaB_z(:,end,:) = [];
    NablaBdata2D.FnablaB_r(end,:,:) = [];
    NablaBdata2D.FnablaB_r(:,end,:) = [];
    NablaBdata2D.FnablaB(end,:,:) = [];
    NablaBdata2D.FnablaB(:,end,:) = [];
    NablaBdata2D.VnablaB_z(end,:,:) = [];
    NablaBdata2D.VnablaB_z(:,end,:) = [];
    NablaBdata2D.VnablaB_r(end,:,:) = [];
    NablaBdata2D.VnablaB_r(:,end,:) = [];
    NablaBdata2D.VnablaB(end,:,:) = [];
    NablaBdata2D.VnablaB(:,end,:) = [];

    save(savename_nablaB,'NablaBdata2D')
end