function [data2D_q,zq,rq] = interp_data(grid2D,data2D,SXRdata)

zq_psi = grid2D.zq;
rq_psi = grid2D.rq;
psi = data2D.psi(:,:,find(data2D.trange==SXRdata.t));
jt = data2D.Jt(:,:,find(data2D.trange==SXRdata.t));
Bz = data2D.Bz(:,:,find(data2D.trange==SXRdata.t));
Br = data2D.Br(:,:,find(data2D.trange==SXRdata.t));

EE = SXRdata.EE;
range = SXRdata.range;

range = range./1000;
zmin1 = range(1);
zmax1 = range(2);
zmin2 = range(3);
zmax2 = range(4);
rmin = range(5);
rmax = range(6);
r_space_SXR = linspace(rmin,rmax,size(EE,1));
z_space_SXR1 = linspace(zmin1,zmax1,size(EE,2));
z_space_SXR2 = linspace(zmin2,zmax2,size(EE,2));
[zq_sxr1,rq_sxr1] = meshgrid(z_space_SXR1,r_space_SXR);
[zq_sxr2,rq_sxr2] = meshgrid(z_space_SXR2,r_space_SXR);

z_space_SXR = linspace(max(zmin1,zmin2),min(zmax1,zmax2),size(EE,2));
[zq_sxr,rq_sxr] = meshgrid(z_space_SXR,r_space_SXR);

for i = 1:4
    if i <= 2
        zq_old = zq_sxr2;
        rq_old = rq_sxr2;
    else
        zq_old = zq_sxr1;
        rq_old = rq_sxr1;
    end
    EE(:,:,i) = interp2(zq_old,rq_old,EE(:,:,i),zq_sxr,rq_sxr);
end


% 異なる2次元グリッド間の座標が重複する範囲を切り出す
rq_min = max(min(rq_psi(:)), min(rq_sxr(:)));
rq_max = min(max(rq_psi(:)), max(rq_sxr(:)));
zq_min = max(min(zq_psi(:)), min(zq_sxr(:)));
zq_max = min(max(zq_psi(:)), max(zq_sxr(:)));

[zq, rq] = meshgrid(linspace(zq_min,zq_max,40), ...
   linspace(rq_min,rq_max,40));

data2D_q.psi = interp2(zq_psi,rq_psi,psi,zq,rq);
data2D_q.jt = interp2(zq_psi,rq_psi,jt,zq,rq);
data2D_q.Bz = interp2(zq_psi,rq_psi,Bz,zq,rq);
data2D_q.Br = interp2(zq_psi,rq_psi,Br,zq,rq);
data2D_q.EE = zeros([size(zq),4]);
for i = 1:4
    data2D_q.EE(:,:,i) = interp2(zq_sxr,rq_sxr,EE(:,:,i),zq,rq);
end

end