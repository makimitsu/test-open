function get_B_reconnection(PCB,pathname)

trange = PCB.trange;
[grid2D,data2D] = process_PCBdata_280ch(PCB,pathname);
Br = data2D.Br;
rq = grid2D.rq;
zq = grid2D.zq;

B_reconnection = zeros(1,numel(trange));
mergingRatio = zeros(1,numel(trange));

for i = 1:numel(trange)
    [magaxis,xpoint] = get_axis_x(grid2D,data2D,i);
    if numel(magaxis.r) == 2
        range = rq>=min(magaxis.r)&rq<=max(magaxis.r)&zq>=min(magaxis.z)&zq<=max(magaxis.z);
        Br_t = Br(:,:,i);
        B_reconnection(1,i) = max(Br_t(range),[],'all');
        mergingRatio(1,i) = xpoint.psi/mean(magaxis.psi);
    else
        B_reconnection(1,i) = NaN;
        mergingRatio(1,i) = NaN;
    end
end

figure;plot(trange,B_reconnection);
figure;plot(trange,mergingRatio);

end