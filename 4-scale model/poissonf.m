% poissonf(750,10,350)
% ���ɷֲ�����

function Pnn=poissonf(d,n,NN)
% d����Ԫ�е�������n����Ԫ�е�����������NN����������������е��������涨Ϊ350

Px_tot=0; % �ۼƸ��ʳ�ֵ
m=d/n;    % ���������е�ƽ������m

Pnn(1,1)=exp(-m); % ����x==0ʱ�ĸ���
Px_tot=Px_tot+Pnn; % �����ۼƸ���

for ak=1:NN % ѭ��NN�Σ�������������1������ʼ��ֱ����������NN��
    Pnn(1,ak+1)=Pnn(1,ak)*m/ak; % ������ak-1����ʱ�ĸ��ʼ���ak����ʱ�ĸ��ʣ�ע�⣺Pnn�е���ű�ʵ�ʵ�������1
    Px_tot=Px_tot+Pnn(1,ak+1); % �����ۼƸ���
end

if Px_tot<0.95 % �������ȫ����Ӻ�С��0.95������ʾ������Ϣ
    msgbox('Possible error in Neyman distribution');
end

for k=1:(NN+1) % ���»��ָ��ʣ�ʹ֮�ֲ���0-1֮��
    Pnn(1,k)=Pnn(1,k)/Px_tot;
end