function result = trigrid_interpor3(z,x,y,yy)

Nterms = length(x)*length(y);
Xtemp = zeros(1,Nterms);
Ytemp = zeros(1,Nterms);
SigTemp = zeros(1,Nterms);

for i = 1:Nterms
    index = mod(i-1,length(y))+1;
    Xtemp(i) = x(floor((i-1)/length(y))+1);
    Ytemp(i) = y(index);
    SigTemp(i) = z(floor((i-1)/length(y))+1,index);
end
Xtemp = Xtemp';
Ytemp = Ytemp';
SigTemp = SigTemp';

F = scatteredInterpolant(Xtemp,Ytemp,SigTemp);
x_return = linspace(min(x),max(x),length(yy));
y_return = y;
[xq,yq] = meshgrid(x_return,y_return);
F.Method = 'linear';
vq = F(xq',yq');

result = struct('x',x_return,'y',y_return,'z',vq);

end