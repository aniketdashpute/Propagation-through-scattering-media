%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code integrates both the steps 1 and 2 to give the final output of the beam
%It iterates over the M screens. 
%My_Step2PassingThroughScreens is the code which runs for every screen
%E1 = Step12Integrator(1,..). Efinal for n=1
%E0 = Step12Integrator(0,..). Efinal for n=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Efinal = My_Step12Integrator(n,h1,map,M,D,L,deltax,factor,Exx,Exy)% Here n is 0 or 1
%deltax is Distance between two sampling points.
N=ceil(D/deltax)*factor+1; % Number of sampling points.
W = D/5;
d = L/M;

x=linspace(-D/2,D/2,N); % Spatial position of the diffuser.
[x, y] = meshgrid(x);

Einitial = ((x+1i*y).^n).*exp(-(x.^2 + y.^2)./(2*(W^2)));
for i=1:N
    for j=1:N
        v = sqrt((i-N/2)^2+(j-N/2)^2);
        % ensure that outside this circle, E is 0
        if (v>N/2)
            Einitial(i,j) = 0;
        end
    end
end
Einitial = Einitial.*map;
% % normalization ??
% Exx = Exx/norm(Einitial,'fro');
% Exy = Exy/norm(Einitial,'fro');
% Einitial = Einitial/norm(Einitial,'fro');

str = sprintf('Initial LG beam n=%d',n);
figure, imagesc(abs(Einitial)), title(str);
%save Einitial

n = mod(n+1,2);
Einitial2 = ((x+1i*y).^n).*exp(-(x.^2 + y.^2)./(2*(W^2)));
for i=1:N
    for j=1:N
        v = sqrt((i-N/2)^2+(j-N/2)^2);
        if (v>N/2)
            Einitial2(i,j) = 0;
        end
    end
end
Einitial2 = Einitial2.*map;
% Einitial2 = Einitial2/norm(Einitial2,'fro');
Ein2 = Einitial2;
str = sprintf('Initial LG beam n=%d',n);
figure, imagesc(abs(Einitial)), title(str);
%save Einitial

%h1 = Step1ScreenParticles(N); %Uncomment to get h1 from here
Ein = Einitial;
for i1=1:M
    S = h1(i1,:,:);
    S = reshape(S,N,N);
    Ein = My_Step2PassingThroughScreens(n,Ein,Ein2,S,M,D,L,deltax,factor,Exx,Exy);
%     keyboard;
%     if (mod(i1,50)==0)
%         display(i1);
%     end
    i1
end

%load Einitial
Efinal = Ein;
str = sprintf('LG beam after passing through %d Screens',M);
% figure, imagesc(abs(Efinal)), title(str);

Efinal = (Efinal);
end