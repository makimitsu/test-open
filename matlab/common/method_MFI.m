function [F] = method_MFI(W,P,n_Vx,n_Vy,d_Vx,d_Vy,plot_analisis)

%‰Šú‰ð‚ð¶¬
[F,ganma] = method_Tikhonov1(W,P,n_Vx,n_Vy,d_Vx,d_Vy,plot_analisis);
filename = ['grad_',num2str(n_Vx),'_',num2str(round(d_Vx,1)),'.mat'];
load(filename ,'g_Vx','g_Vy')
disp('Saved g_Vx and g_Vy is loaded.')

j = 1;
res = 1;
if plot_analisis
    figure('Position',[300 150 550 500])
end

while res > 0.1 && j<50
    disp(j)
    disp(res)
    if plot_analisis
        plot(j,res, 'r+')
        hold on
    end
    buf = F;
    M = zeros(n_Vx*n_Vy);
    for k = 1:n_Vx*n_Vy
        if F(k) > 0
            M(k,k)=1/F(k);
        end
    end
    max_M = max(max(M));
    for k = 1:n_Vx*n_Vy
        if F(k) <= 0
            M(k,k)=max_M;
        end
    end
    penalty = g_Vx.'*M*g_Vx + g_Vy.'*M*g_Vy;
    F = inv(W.'*W + ganma * penalty)*W.'*P;
    res = norm(F - buf);
    j = j+1;
end
if plot_analisis
    xlabel('Iteration count','FontSize',15);
    ylabel('Residual','FontSize',15);
    hold off
end
F = method_goto0(F);
end
