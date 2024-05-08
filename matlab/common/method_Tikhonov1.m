function [F,ganma] = method_Tikhonov1(W,P,n_Vx,n_Vy,d_Vx,d_Vy,plot_analisis)
max_i = 100;
p = zeros(2,max_i);
curv = zeros(2,max_i);%1s–ÚƒKƒ“ƒ}A2s–Ú‹È—¦
filename = ['grad_',num2str(n_Vx),'_',num2str(round(d_Vx,1)),'.mat'];
if exist(filename ,'file')
    load(filename ,'grad')
    disp('Saved grad is loaded.')
else
    g_Vx = zeros(n_Vx*n_Vy);
    g_Vy = zeros(n_Vx*n_Vy);
    for i = 1:n_Vx*n_Vy
        x1 = 1 + idivide(i-1, int8(n_Vy), 'floor');
        y1 = int8(mod(i,n_Vy));
        if x1 ~= 1 && x1 ~= n_Vx
            for j = 1:n_Vx*n_Vy
                x2 = 1 + idivide(j-1, int8(n_Vy), 'floor');
                y2 = int8(mod(j,n_Vy));
                if y1 == y2
                    if x2 == x1
                        g_Vx(i,j) = -1/d_Vx;
                    elseif x2 == x1 + 1
                        g_Vx(i,j) = 1/d_Vx;
                    end
                end
            end
        end
        if y1 ~= 1 && y1 ~= n_Vy
            for j = 1:n_Vx*n_Vy
                x2 = 1 + idivide(j-1, int8(n_Vy), 'floor');
                y2 = int8(mod(j,n_Vy));
                if x1 == x2
                    if y2 == y1
                        g_Vy(i,j) = -1/d_Vy;
                    elseif y2 == y1 + 1
                        g_Vy(i,j) = 1/d_Vy;
                    end
                end
            end
        end
    end
    grad = g_Vx.'*g_Vx + g_Vy.'*g_Vy;
    disp('grad is generated.')
    %grad‚ð•Û‘¶(ŽŸ‰ñˆÈ~ŒvŽZŽžŠÔß–ñ‚Ì‚½‚ß)
    save(filename,'grad','g_Vx','g_Vy')
end
if plot_analisis
    figure('Position',[300 150 550 500])
end
for i = 1:max_i
    disp(i)
    ganma = 1.1^(-i+10);
    curv(1,i) = ganma;
    F = inv(W.'*W + ganma * grad)*W.'*P;
    p(1,i) = log10(norm(W * F - P));
    p(2,i) = log10(norm(F));
    if plot_analisis
        plot(p(1,i),p(2,i), 'r+')
        hold on
    end
    if i>1 && i<max_i
        curv(2,i) = 2*det([p(:,i+1)-p(:,i), p(:,i-1)-p(:,i)])/...
            (norm(p(:,i+1)-p(:,i))*norm(p(:,i-1)-p(:,i))*norm(p(:,i+1)-p(:,i-1)));
    end
end
if plot_analisis
    title('L-curve','FontSize',20);
    xlabel('log10(WF-P)','FontSize',15);
    ylabel('log10(F)','FontSize',15);
    hold off
    figure('Position',[500 150 550 500])
    title('Curvature of L-curve','FontSize',20);
    plot(p(1,:),curv(2,:))
end
[~,I] = max(curv.');
ganma = curv(1,I(2));%Å“K‚ÈƒKƒ“ƒ}‚ð—^‚¦‚é
F = inv(W.'*W + ganma * grad)*W.'*P;
F = method_goto0(F);
end