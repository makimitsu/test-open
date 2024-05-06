function result = trigrid_interpor_for_r_lambda2(z,x,y,lambda_out,yy)

Nterms = length(x)*length(y);
Xtemp = zeros(1,Nterms);
Ytemp = zeros(1,Nterms);
SigTemp = zeros(1,Nterms);

for i = 1:Nterms
    index = mod(i-1,length(y))+1;
    index2 = floor((i-1)/length(y))+1;
    Xtemp(i) = x(index2,index);
    Ytemp(i) = y(index);
    SigTemp(i) = z(index2,index);
end
Xtemp = Xtemp';
Ytemp = Ytemp';
SigTemp = SigTemp';

F = scatteredInterpolant(Xtemp,Ytemp,SigTemp);
% x_interp = linspace(min(x,[],"all"),max(x,[],"all"),length(lambda_out));
% y_interp = linspace(min(y,[],"all"),max(y,[],"all"),length(yy));
% [xq,yq] = meshgrid(x_interp,y_interp);
[xq,yq] = meshgrid(lambda_out,yy);
F.Method = 'linear';
vq = F(xq,yq);

% figure(2)
% subplot(1,2,1)
% contourf(z');
% subplot(1,2,2)
% contourf(xq,yq,vq);
% pause(0.5)

result = struct('xq',xq,'yq',yq,'zq',vq);
end