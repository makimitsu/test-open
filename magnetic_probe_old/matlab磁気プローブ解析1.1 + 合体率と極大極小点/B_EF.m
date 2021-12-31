function [Bz_EF,Br_EF] = B_EF(z_EF,r_EF,i_EF,n_EF,probe_mesh_r,probe_mesh_z)
% input:
%   (z,r) position of center of EF coils
%   EF coil current and turns
%   2d mesh of (z,r) position of B probes
% output:
%   2d array of Bz and Br due to EF coil,
%       and it corresponds to the input position mesh
% a-helmholtz-coil-like configuration is assumed
% algorithm comes from:
%   https://tiggerntatie.github.io/emagnet/offaxis/iloopoffaxis.htm
%   https://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-UR-17-27619

mu0   = 4*pi*10^(-7);
alpha = probe_mesh_r/r_EF;
beta1 = (probe_mesh_z+z_EF)/r_EF;
beta2 = (probe_mesh_z-z_EF)/r_EF;
Q1 = ((1+alpha).^2+beta1.^2);
Q2 = ((1+alpha).^2+beta2.^2);
k1 = (4*alpha./Q1).^0.5;
k2 = (4*alpha./Q2).^0.5;
m1 = k1.^2;
m2 = k2.^2;
B0 = i_EF*n_EF*mu0/(2*r_EF);

Bz1 = B0./(pi*Q1.^0.5).*(ellipticE(m1).*(1-alpha.^2-beta1.^2)./(Q1-4*alpha)+ellipticK(m1));
Bz2 = B0./(pi*Q2.^0.5).*(ellipticE(m2).*(1-alpha.^2-beta2.^2)./(Q2-4*alpha)+ellipticK(m2));
Bz_EF = Bz1 + Bz2;

gamma1 = (probe_mesh_z+z_EF)./probe_mesh_r;
gamma2 = (probe_mesh_z-z_EF)./probe_mesh_r;
Br1 = B0*gamma1./(pi*Q1.^0.5).*(ellipticE(m1).*(1+alpha.^2+beta1.^2)./(Q1-4*alpha)-ellipticK(m1));
Br2 = B0*gamma2./(pi*Q2.^0.5).*(ellipticE(m2).*(1+alpha.^2+beta2.^2)./(Q2-4*alpha)-ellipticK(m2));
Br_EF = Br1 + Br2;
end