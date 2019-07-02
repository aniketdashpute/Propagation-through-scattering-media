function [] = main_particles
clear all
for k=1:20
    G = 100*k;
    [E0,E1,E,Ef0,Ef1,Ef,h] = main_My_Mie_particles(G);
    s1 = 'Mie_data_';
    s2 = num2str(G);
    str = strcat(s1,s2);
    save(str);
%     keyboard;
end
close all
end

