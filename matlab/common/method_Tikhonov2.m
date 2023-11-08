function [F] = method_Tikhonov2(W,P,n_Vx,n_Vy,d_Vx,d_Vy,plot_analisis)
max_i = 10;
p = zeros(2,max_i);
curv = zeros(2,max_i);%1s–ÚƒKƒ“ƒ}A2s–Ú‹È—¦
filename = ['gradgrad_',num2str(n_Vx),'_',num2str(round(d_Vx,1)),'.mat'];
if exist(filename ,'file')
    load(filename ,'gradgrad')
    disp('Saved gradgrad is loaded.')
else
    gradgradVx = zeros(n_Vx*n_Vy);
    gradgradVy = zeros(n_Vx*n_Vy);
    for i = 1:n_Vx*n_Vy
        x1 = 1 + idivide(i-1, int8(n_Vy), 'floor');
        y1 = int8(mod(i,n_Vy));
        if x1 ~= 1 && x1 ~= n_Vx
            for j = 1:n_Vx*n_Vy
                x2 = 1 + idivide(j-1, int8(n_Vy), 'floor');
                y2 = int8(mod(j,n_Vy));
                if y1 == y2
                    if x2 == x1 - 1
                        gradgradVx(i,j) = 1/d_Vx^2;
                    elseif x2 == x1
                        gradgradVx(i,j) = -2/d_Vx^2;
                    elseif x2 == x1 + 1
                        gradgradVx(i,j) = 1/d_Vx^2;
                    end
                end
            end
        end
        if y1 ~= 1 && y1 ~= n_Vy
            for j = 1:n_Vx*n_Vy
                x2 = 1 + idivide(j-1, int8(n_Vy), 'floor');
                y2 = int8(mod(j,n_Vy));
                if x1 == x2
                    if y2 == y1 - 1
                        gradgradVy(i,j) = 1/d_Vy^2;
                    elseif y2 == y1
                        gradgradVy(i,j) = -2/d_Vy^2;
                    elseif y2 == y1 + 1
                        gradgradVy(i,j) = 1/d_Vy^2;
                    end
                end
            end
        end
    end
    gradgrad = gradgradVx.'*gradgradVx + gradgradVy.'*gradgradVy;
    disp('gradgrad is generated.')
    %grad‚ð•Û‘¶(ŽŸ‰ñˆÈ~ŒvŽZŽžŠÔß–ñ‚Ì‚½‚ß)
    save(filename,'gradgrad')
end
if plot_analisis
    figure('Position',[300 150 550 500])
end
for i = 1:max_i
    disp(i)
    ganma = 1.1^(-i+10);
    curv(1,i) = ganma;
    F = inv(W.'*W + ganma * gradgrad)*W.'*P;
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
    ylabel('log10(Þ^2F)','FontSize',15);
    hold off
    figure('Position',[500 150 550 500])
    title('Curvature of L-curve','FontSize',20);
    plot(p(1,:),curv(2,:))
end
[~,I] = max(curv.');
ganma = curv(1,I(2));%Å“K‚ÈƒKƒ“ƒ}‚ð—^‚¦‚é
F = inv(W.'*W + ganma * gradgrad)*W.'*P;
F = method_goto0(F);
end