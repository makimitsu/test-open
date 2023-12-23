clear;

n=90;
z_0=0;
r_0=-0.3;
z=linspace(-1,1,n);
r=linspace(-1,1,n);

[z_space,r_space] = meshgrid(z,r); %rが横、zが縦の座標系（左上最小）
EE = exp(-(r_space.^2+z_space.^2));
EE = flipud(EE);

figure;imagesc(EE);colorbar;

E = reshape(EE,[],1);

[grad_r,grad_z] = grad_matrix(n);
Er = grad_r*E;
Ez = grad_z*E;

EEr = reshape(Er,n,n);
EEz = reshape(Ez,n,n);

figure;
subplot(1,2,1);imagesc(EEr);title('r gradient');
subplot(1,2,2);imagesc(EEz);title('z gradient');