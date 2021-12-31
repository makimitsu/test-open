function psi = get_psi(B_z,r_probe,time)
% input:
%   3d array of double: B_z (r,z,t)
%   1d array of double: r_probe, locations of probes along r
%   integer: time of interest
% output:
%   2d array of double: psi (r,z), the poloidal magnetic flux on rz plane

B_z_t = B_z(:,:,time);
[r_count,z_count] = size(B_z_t);
psi = zeros(r_count,z_count);

r_step = 0.001;
r_interp = r_probe(1):-r_step:r_probe(end);
epsilon = 10^(-10);
r_index = r_probe;
for i = 1:r_count
    r_index(i) = find(abs(r_interp - r_probe(i)) < epsilon);
end

B_z_interp = zeros(length(r_interp),z_count);

for j = 1:z_count
    psi(r_count,j) = 2*pi*r_probe(r_count)*B_z_t(r_count,j)*r_probe(r_count);
    B_z_interp(:,j) = interp1(r_probe,B_z_t(:,j),r_interp,'linear');
end

r_interp = rot90(r_interp,3);
B_z_times_r = B_z_interp.*r_interp;

for i = r_count-1:-1:1
    for j = 1:z_count
        psi(i,j) = psi(i+1,j)+2*pi*ones(1,r_index(i+1)-r_index(i))*r_step*B_z_times_r(r_index(i):r_index(i+1)-1,j);
    end
end

end

