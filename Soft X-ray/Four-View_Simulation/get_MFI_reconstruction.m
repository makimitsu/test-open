function g = get_MFI_reconstruction(f,T,N_grid)
% f: vector of projected data (input)
% g: vector of local distribution (output)
% T: geometric weight matrix (known)
% N_grid: grid number of local distribution

sigma = std(f);
T_ = T./sigma;
f_ = f./sigma;

[grad_r,grad_z] = grad_matrix(N_grid);
H = grad_r.'*grad_r + grad_z.'*grad_z; %H0

W = zeros(size(grad_r));
n = 3; % max number of repetition
for i = 1:n
    % g = inv(T_.'*T_+a*H)*(T_.'*f_);
    g = (T_.'*T_+a*H)\(T_.'*f_);
    for j = 1:numel(diag(W))
        if g(j) > 0
            W(j,j) = 1./g(j);
        else
            W(j,j) = -1;
        end
        W(W==-1) = max(W);
    end
    H = grad_r.'*W*grad_r + grad_z.'*W*grad_z; %H_(n)
end

end

function [grad_r,grad_z] = grad_matrix(N_grid)

grad_r = zeros(N_grid^2);
grad_z = grad_r;

A = zeros(N_grid);
for i = 1:N_grid
    if 1<i && i<N_grid
        A(i,i) = 0;
        A(i,i-1) = 1/2;
        A(i,i+1) = -1/2;
    elseif i == 1
        A(i,i) = 1/2;
        A(i,i+1) = -1/2;
    else
        A(i,i) = -1/2;
        A(i,i-1) = 1/2;
    end
end

for i = 1:N_grid
    grad_r(N_grid*(i-1)+1:N_grid*i,N_grid*(i-1)+1:N_grid*i) = A;
end

for i = 1:N_grid^2
    if i <= N_grid
        grad_z(i,i) = -1/2;
        grad_z(i,i+N_grid) = 1/2;
    elseif i <= N_grid*(N_grid-1)
        grad_z(i,i) = 0;
        grad_z(i,i-N_grid) = -1/2;
        grad_z(i,i+N_grid) = 1/2;
    else
        grad_z(i,i) = 1/2;
        grad_z(i,i-N_grid) = -1/2;
    end
end

end