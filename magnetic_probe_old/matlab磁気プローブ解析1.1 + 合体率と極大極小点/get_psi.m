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

for j = 1:z_count
    psi(r_count,j) = 2*pi*r_probe(r_count)*B_z_t(r_count,j)*r_probe(r_count);
end

for i = r_count-1:-1:1
    for j = 1:z_count
        psi(i,j) = psi(i+1,j)+2*pi*r_probe(i)*(B_z_t(i,j)+B_z_t(i+1,j))/2*(r_probe(i)-r_probe(i+1));
    end
end

end

