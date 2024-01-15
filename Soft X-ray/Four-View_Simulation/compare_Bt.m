function [] = compare_Bt(PCB,pathname,time)

% [grid2D,data2D] = process_PCBdata_280ch(PCB,pathname);
trange = PCB.trange;
PCB_ref = PCB;
PCB_ref.shot = PCB.tfshot;
PCB_ref.tfshot = [0,0];
[grid2D,data2D_ref] = process_PCBdata_280ch(PCB_ref,pathname);
Bt_ref = data2D_ref.Bt;

[I_TF,x,aquisition_rate] = get_TF_current(PCB,pathname);
m0 = 4*pi*10^(-7);
timing = x/aquisition_rate==time;
B_guide = m0*I_TF(timing)*1e3*12./(2*pi()*grid2D.rq);


figure('Position', [0 0 1500 1500],'visible','on');
subplot(1,2,1);
contourf(grid2D.zq(1,:),grid2D.rq(:,1),Bt_ref(:,:,trange==time),40,'LineStyle','none')
subplot(1,2,2);
contourf(grid2D.zq(1,:),grid2D.rq(:,1),B_guide,40,'LineStyle','none')

figure;imagesc(Bt_ref(:,:,trange==time)./B_guide);
clim([0 1]);colorbar;

end