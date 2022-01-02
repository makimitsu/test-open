function [Et] = get_Et(B_z,r_probe)

psi_all = zeros(size(B_z));

for i = 1:1024
    psi_all(:,:,i) = get_psi(B_z,r_probe,i);
end

for i = 1:length(r_probe)
    Et = -1/(2*pi*r_probe(i))*diff(psi_all,1,3);
end
end

%     psi_before = get_psi(B_z,r_probe,time-1);
%     psi_now = get_psi(B_z,r_probe,time);
%     E_t = zeros(size(psi_now));
%     for i = 1:length(r_probe)
%         E_t(i,:) = -1/(2*pi*r_probe(i)) * (psi_now(i,:) - psi_before(i,:));
%     end
% end