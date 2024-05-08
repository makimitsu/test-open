function [F] = method_SIRT(W,P,plot_analisis)
[~,n_cols] = size(W);
rho = 10^(-3);%ŠÉ˜aƒpƒ‰ƒ[ƒ^
j = 1; 
res = 1;
F = zeros(n_cols,1);
if plot_analisis
    figure('Position',[300 150 550 500])
end
while res > 1*10^-2 && j<100
    buf = F;
    F = F + rho*W.'*(P-W*F);
    F = method_goto0(F);
    res = norm(F - buf);
    if plot_analisis
        plot(j,res, 'r+')
        hold on
    end
    j = j+1;
end
if plot_analisis
    xlabel('Iteration count','FontSize',15);
    ylabel('Residual','FontSize',15);
    hold off
end
