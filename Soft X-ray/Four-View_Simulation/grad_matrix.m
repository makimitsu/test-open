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