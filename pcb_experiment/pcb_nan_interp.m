function vq = pcb_nan_interp(rpos_bt, zpos_bt, grid2D, bt, ok_bt, t)

r_bt_probe = rpos_bt(11:20);
z_bt_probe = zpos_bt(2:10:92);
rpos_bt(4:10) = rpos_bt(14:20);
rpos_bt(21) = rpos_bt(11);
rpos_bt(94:100) = rpos_bt(14:20);
zpos_bt(4:10) = -0.17*ones(size(zpos_bt(4:10)));
zpos_bt(21) = -0.085;
zpos_bt(94:100) = 0.17*ones(size(zpos_bt(94:100)));

vq = zeros(10,10);

for i = 1:100
    ir_bt = find(r_bt_probe == rpos_bt(i));
    iz_bt = find(z_bt_probe == zpos_bt(i));
    if (ok_bt(i))
        vq(ir_bt,iz_bt) = bt(t,i) - mean(bt(10:40,i));
    else
        vq(ir_bt,iz_bt) = NaN;
    end
end

vq = inpaint_nans(vq,0);
vq = griddata(z_bt_probe,r_bt_probe,vq,grid2D.zq,grid2D.rq);

end