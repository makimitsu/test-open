x = linspace(-100,100,201);
v1 = exp(-x.^2/1000);%���u�֐�
v2 = v1(80:120);
% q = 2*exp(-x.^2/100)+exp(-(x-10).^2/1000);%�^�̃X�y�N�g��
q = 2*exp(-(x+30).^2/100);%�^�̃X�y�N�g��
u1 = conv(q,v2);%�ϑ��X�y�N�g��
u2 = conv(q,v1,'same');%�ϑ��X�y�N�g��
[w,r] = deconv(u1,v2); 
u3 = conv(w,v1,'same');
figure
% plot(x,v1,'k',x,q,'r',x,u2/10,'g')
% legend('���u�֐�','�^�̃X�y�N�g��','�ϑ��X�y�N�g��')
% plot(x,v1,'k',x,q,'r',x,w,'b')
% legend('���u�֐�','�^�̃X�y�N�g��','�����X�y�N�g��')
plot(x,v1,'k',x,q+0.1,'r',x,u2/16+0.1,'g',x,w,'b',x,u3/16,'m')
legend('���u�֐�','�^�̃X�y�N�g��','�ϑ��X�y�N�g��','�����X�y�N�g��','�ϑ��X�y�N�g��(�Č�)')
