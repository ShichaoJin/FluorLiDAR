% neyman(3000,40,20,350);
% Ŧ��A�ֲ�
function Px=neyman(d,n,m2,NN)
% d����Ԫ�е�������n����Ԫ�е�����������m2��һ��Ⱥ���ƽ��������NN����������������е��������涨Ϊ350

m1=d/(n*m2); % ����m1��һ������ƽ���м���Ⱥ��

Px(1,1)=exp(-m1*(1-exp(-m2))); % ��������0�����ĸ���,������Px�����е������1
Px_tot=Px(1,1); % �ۼƸ���
for k=2:(NN+1) % ��������1������NN����ʱ�ĸ���
    Px(1,k)=0; % ����ʼֵ
    if (k-1)>d/n && Px(1,k-1)<0.00001 % Ϊ�ӿ��ٶ�,�������������ص�ƽ��ֵ�Ҹ���С��0.00001�Ĳ�����
		Px(1,k)=0; 
    else
        for t=1:(k-1)
            PXX=m1*m2*exp(-m2)/(k-1)*Px(1,k-t);
            for i=2:t
                PXX=PXX*m2/(i-1);
            end
            Px(1,k)=Px(1,k)+PXX;
        end
    end
    Px_tot=Px_tot+Px(1,k); %  making sure that the sum of Px is unity
end

if Px_tot<0.98 % �������ȫ����Ӻ�С��0.98������ʾ������Ϣ
    msgbox('Possible error in Neyman distribution');
end

for k=1:(NN+1) % ���»��ָ��ʣ�ʹ֮�ֲ���0-1֮��
    Px(1,k)=Px(1,k)/Px_tot;
end