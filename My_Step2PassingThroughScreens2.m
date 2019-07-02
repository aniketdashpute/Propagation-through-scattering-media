%W = D/5;
%E = ((x+y*1i)^n)*exp(-(x^2+y^2)/(2*W^2))*exp(1i*kz*z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code propagates the beam after passing through the screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Eout = My_Step2PassingThroughScreens2(n,Ein,S,M,D,L)
%%----------------- Old Code Values---------------------%%
deltax = 0.2;
N=ceil(D/deltax)+1; % Number of sampling points.
d = L/M;
n0=1.025; % Relative refractive index of the particle to the medium.
n1=1.5; % Refractive index of the medium.
mean=0.2; % Particle diameter.

lamda = 0.5;
kz = 2*pi/lamda;
lamda1 = lamda;

k0 = kz; % vacuum wave number
a=deltax/2; % sphere radius

fx=1/deltax*(-0.5:1/(N-1):0.5); % Spatial frequencies.
[fx, fy] = meshgrid(fx);
x=linspace(-D/2,D/2,N); % Spatial position of the diffuser.
[x, y] = meshgrid(x);

%%-------------- Calculations -------------------------%%
[m1,n1] = size(S); %Size of the matrix which has scatterers
% [m2,n2] = size(S); %Size of the matrix on the receving end. Here it is kept the same

% z=1; %Will have to set it
z=d;
Ex = Ein;
Ex_new = zeros(size(S)); %Size of the matrix on the receving end
    for i=1:m1
        for j=1:n1
            if (S(i,j)==1)
                x0 = x(i,j);
                y0 = y(i,j);
                R = ((x0-x).^2+(y0-y).^2+(z)^2).^0.5;
                CTh = z./R;
%                  % Make cosTh matrix for a scattering pixel
%                 CTh = zeros(size(I2));
%                 for a=1:m2
%                     for b=1:n2
%                         x0 = x(i,j);
%                         y0 = y(i,j);
%                         xd = x(a,b);
%                         yd = y(a,b);
%                         R = ((x0-xd)^2+(y0-yd)^2+(z)^2)^0.5;
%                         cosTh = z/R;
%                         CTh(a,b) = cosTh;
%                     end
%                 end
                m=n0;
                % size parameter x=k0*a, and u=cos(scattering angle),
                % where k0=vacuum wave number, a=sphere radius;
                xx=k0*a;
                u=CTh;
                Ex_new_temp = f1(n,CTh,m,xx,u,z,kz).*Ex;
                Ex_new = Ex_new + Ex_new_temp;
            end
        end
        i
    end    
Escat = Ex_new;

%-----------------------Add direct beam------------------%
H2=fft2(fftshift(Ein)); %H2 is the Fourier transform of the input field.
P2=1.*(sqrt(fx.^2+fy.^2)<=n1/lamda1)-1i.*(sqrt(fx.^2+fy.^2)>n1/lamda1);% This step assigns a value of 1 for ?fx2 + fy2 <= n1/lamda1 and of i for ?fx2 + fy2 >n 1/lamda1.
P2=exp(-1i*2*pi*d.*P2.*sqrt(abs((n1/lamda1)^2-fx.^2-fy.^2)));%Transfer function of free space propagation, see Eq. (13).
P2=fftshift(P2); %equation
G2=H2.*P2; % Direct beam calculation executed at the Fourier domain.
Edirect=ifftshift(ifft2(G2)); % Retrieving direct beam with inverse Fourier transform.
%--------------------------------------------------------%

Eout = Escat+Edirect;
end



% function [Efx] = f1(n,cosTh,m,x,u,z,kz)
% R = z./cosTh;
% [m1,n1] = size(u);
% S12 = zeros(size(u));
% for i=1:m1
%     for j=1:n1
%         ud = u(i,j);
%         S = Mie_S12(m,x,ud);
%         S12(i,j) = S(n+1);
%     end
% end
% % Efx = exp(-1i*k*z/cosTh)*S2;
% Efx = (exp(-1i*kz*R+1i*kz*z)./(1i*kz*R)).*S12;
% end

function [Efx] = f1(n,cosTh,m,x,u,z,kz)
R = z./cosTh;
% u=1;
S = Mie_S12(m,x,u,n);
S12 = S(n+1);
% Efx = exp(-1i*k*z/cosTh)*S2;
Efx = (exp(-1i*kz*R+1i*kz*z)./(1i*kz*R)).*S12;
end



















% H1=fft2(fftshift(Ein.*S)); %H1 is the Fourier transform of the input field times the 2D particle position map.
% H2=fft2(fftshift(Ein)); %H2 is the Fourier transform of the input field.
% 
% P2=1.*(sqrt(fx.^2+fy.^2)<=n1/lamda1)-1i.*(sqrt(fx.^2+fy.^2)>n1/lamda1);% This step assigns a value of 1 for ?fx2 + fy2 <= n1/lamda1 and of i for ?fx2 + fy2 >n 1/lamda1.
% P2=exp(-1i*2*pi*d.*P2.*sqrt(abs((n1/lamda1)^2-fx.^2-fy.^2)));%Transfer function of free space propagation, see Eq. (13).
% P2=fftshift(P2); %equation.
% 
% P1=n1^2*(2*pi/lamda1)^2*(mean/2)^3*(n^2-1)/(n^2+2)*(y.^2+d^2).*exp(1i*2*pi*n1/lamda1*sqrt(x.^2+y.^2+d^2))./sqrt(x.^2+y.^2+d^2).^3; %Scattering equation from Eq. (11).
% P1=fft2(fftshift(P1)); % Fourier transform of the scattering formula.
% G1=H1.*P1; % Scattering calculation executed at the Fourier domain.
% G2=H2.*P2; % Direct beam calculation executed at the Fourier domain.
% g1=ifftshift(ifft2(G1)); % Retrieving scattering field with inverse Fourier transform.
% g3=ifftshift(ifft2(G2)); % Retrieving direct beam with inverse Fourier transform.
% Eprop=g3+g1; % The output field is the sum of the scattering field and the direct beam.
% % Eprop = Eprop.*exp(1i*kz*d);
% Eout = Eprop;
