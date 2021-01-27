function [SNR1,SNR2,SNR3,x] = main_SNR_particles(r)
for k=1:5
    disp(k);
    M = k;
    [s1,s2,s3] = main_SNR(M);
    close all;
    % calculate SNR for the given radius
    SNR1(k) = s1(r);
    SNR2(k) = s2(r);
    SNR3(k) = s3(r);
    x(k)=k;
end

s = 'Radius = ';
rr = num2str(r);
str = strcat(s,rr);
disp(str);

% Plots
figure, plot(x,SNR1,'Color','b'), title(str);
hold on
plot(x,SNR2,'Color','r');
hold off
hold on
plot(x,SNR3,'Color','g');
hold off
end