function g = get_MFI_reconstruction(f,T)
% f: vector of projected data (input)
% g: vector of local distribution (output)
% T: geometric weight matrix (known)
% N_grid: grid number of local distribution

f = f.'; %fを縦ベクトルに変換
% sigma = 5; % error of experimental (projected) data
e_min = 1; % minimum error of experimental (projected) data

snr = 0.1; %expected signal-to-noise ratio

N_grid = sqrt(size(T,2));
[grad_r,grad_z] = grad_matrix(N_grid);
H = grad_r.'*grad_r + grad_z.'*grad_z; %H0

% a = 1;

W = eye(size(grad_r));
diag_idx = find(W);
n = 5; % max repetition number of outer loop
m = 5; % max repetition number of inner loop
a = ones(1,m);
for i = 1:n
    if i > 1
        g_old = g;
    end
    for j = 1:m
        sigma = f*snr;
        sigma(sigma<e_min) = e_min;
        g = (T.'*T+a(j)*H)\(T.'*f);
        Chi_square = mean(((f-T*g)./sigma).^2);
        % Chi_square = mean(((f-T*g)/sigma).^2);
        disp(strcat('iteration : ',num2str(i),'-',num2str(j)));
        disp(strcat('alpha : ',num2str(a(j))));
        disp(strcat('Chi-square : ',num2str(Chi_square)));
        % 再構成結果の誤差が計測誤差程度になれば終了
        if abs(Chi_square-1) <= 0.3
            % disp(a(j));
            break
        elseif j<m 
            if Chi_square < 1
                if j > 1 && a(j) < a(j-1)
                    a(j+1) = (a(j)+a(j-1))/2;
                else
                    a(j+1) = a(j)*10;
                end
            else
                if j > 1 && a(j) > a(j-1)
                    a(j+1) = (a(j)+a(j-1))/2;
                else
                    a(j+1) = a(j)/10;
                end
            end
        end
    end
    if i < n
        g_sample = g;
        g_sample(g_sample<=0) = -1;
        W(diag_idx) = 1./g_sample;
        W(W<0) = max(W(diag_idx));
        H = grad_r.'*W*grad_r + grad_z.'*W*grad_z; %H_(n)
    end
    if i > 1
        figure;
        plot(g_old);
        hold on
        plot(g);
        legend({'g(n-1)','g(n)'});
    end
end

% g(g<0) = 0;

end

% function [grad_r,grad_z] = grad_matrix(N_grid)
% 
% grad_r = zeros(N_grid^2);
% grad_z = grad_r;
% 
% A = zeros(N_grid);
% for i = 1:N_grid
%     if 1<i && i<N_grid
%         A(i,i) = 0;
%         A(i,i-1) = 1/2;
%         A(i,i+1) = -1/2;
%     elseif i == 1
%         A(i,i) = 1/2;
%         A(i,i+1) = -1/2;
%     else
%         A(i,i) = -1/2;
%         A(i,i-1) = 1/2;
%     end
% end
% 
% for i = 1:N_grid
%     grad_r(N_grid*(i-1)+1:N_grid*i,N_grid*(i-1)+1:N_grid*i) = A;
% end
% 
% for i = 1:N_grid^2
%     if i <= N_grid
%         grad_z(i,i) = -1/2;
%         grad_z(i,i+N_grid) = 1/2;
%     elseif i <= N_grid*(N_grid-1)
%         grad_z(i,i) = 0;
%         grad_z(i,i-N_grid) = -1/2;
%         grad_z(i,i+N_grid) = 1/2;
%     else
%         grad_z(i,i) = 1/2;
%         grad_z(i,i-N_grid) = -1/2;
%     end
% end
% 
% end