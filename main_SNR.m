function [SNR1,SNR2,SNR3,x] = main_SNR(screens)
% s1 = 'sd_Mie_data_';
s1 = 'phScr_Mie_data_';
s2 = num2str(screens);
str = strcat(s1,s2);
load(str);
s1 = 'Ef';
s2 = 'Ef0';
s3 = 'Ef1';
% str1 = strcat(s1,str);
% str2 = strcat(s2,str);
% str3 = strcat(s3,str);

I1 = abs(eval(s1));
I2 = abs(eval(s2));
I3 = abs(eval(s3));
[n,m] = size(I1);

[max_num,max_idx] = max(I1(:));
[X,Y]=ind2sub(size(I1),max_idx);
[ll,bb] = size(I1);
X = floor(ll/2);
Y = floor(bb/2);
% disp(X);
% disp(Y);

D = 1;
HM = 40;
for p=1:HM
%     C = ( (x.^2)+(y.^2) -(D/2).^2 );
%     C1 = 1.*(C<=0)+0.*(C>0);
%     I1 = I1.*C1;
%     I2 = I2.*C2;
%     I3 = I3.*C3;
    k=0;
    data1(1)=0;
    data2(1)=0;
    data3(1)=0;
    for i=1:n
        for j=1:m
            d = ((X-i)^2+(Y-j)^2)^0.5;
            if (d<=D)
                k=k+1;
                data1(k) = I1(i,j);
                data2(k) = I2(i,j);
                data3(k) = I3(i,j);            
            end
        end
    end
%     m1 = mean2(I1);
    m1 = mean2(data1);
    sd1 = std2(data1);
%     m2 = mean2(I2);
    m2 = mean2(data2);
    sd2 = std2(data2);
%     m3 = mean2(I3);
    m3 = mean2(data3);
    sd3 = std2(data3);
    SNR1(p) = m1/sd1;
    SNR2(p) = m2/sd2;
    SNR3(p) = m3/sd3;
    x(p) = D;        
    D = D+1;
end

% s = 'Tapes = ';
% ss = strcat(s,str);
figure, plot(x,SNR1,'Color','b');
hold on
plot(x,SNR2,'Color','r');
hold off
hold on
plot(x,SNR3,'Color','g');
hold off

SNR1 = reshape(SNR1,[HM,1]);
SNR2 = reshape(SNR2,[HM,1]);
SNR3 = reshape(SNR3,[HM,1]);
x = reshape(x,[HM,1]);

end