%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code propagates the beam after passing through the screen
%In the end it adds both the components of beam, direct and scattered
%Direct beam is obtained after propagation
%Scattered beam is obtained using Mie scattering

% Ein is input to the screen in consideration
% S is the distribution of particles in this screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Eout = My_Step2PassingThroughScreens(n,Ein,Ein2,S,M,D,L,deltax,factor,Exx,Exy)
%%----------------- Old Code Values---------------------%%
N=ceil(D/deltax)*factor+1; % Number of sampling points.
d = L/M;
n0=1.025; % Relative refractive index of the particle to the medium.
n1=1.5; % Refractive index of the medium.
% mean=0.2; % Particle diameter.

lamda = 0.5;
kz = 2*pi/lamda;
lamda1 = lamda;

k0 = kz; % vacuum wave number
a=deltax/2; % sphere radius

fx=1/deltax*(-0.5:1/(N-1):0.5)*factor; % Spatial frequencies.
[fx, fy] = meshgrid(fx);
x=linspace(-D/2,D/2,N); % Spatial position of the diffuser.
[x, y] = meshgrid(x);

%%-------------- Calculations -------------------------%%
% [m1,n1] = size(S); %Size of the matrix which has scatterers
% [m2,n2] = size(S); %Size of the matrix on the receving end. Here it is kept the same

z=d;

%%---------Response to be used for convolution----------%%
    x0 = 0;
    y0 = 0;
    R = ((x0-x).^2+(y0-y).^2+(z)^2).^0.5;
    CTh = z./R; %Matrix with all cos values
    m=n0;
    % for complex refractive index m=m'+im", 
    % size parameter xx=k0*a, and u=cos(scattering angle),
    % where k0=vacuum wave number, a=sphere radius;
    xx=k0*a;
    u=CTh;
	
	% f1 is a function defined below
%     Ex_response = f1(n,CTh,m,xx,u,z,kz);
%     figure, imagesc(abs(Ex_response));

%---------------Convolution----------------------------%
%------ Scattering part due to the wave passing -------%
Ex_response = Exx;
F_Ein=fftshift(fft2(ifftshift(S.*Ein))); %Fourier transform of the input field.
figure, imagesc(abs(Exx)), title('Exx');

F_response = fftshift(fft2(ifftshift(Ex_response))); %Fourier transform of the response
% figure, imagesc(abs(F_response));

F_Escat = F_Ein.*F_response; %Direct scattered field calculation
Escat=fftshift(ifft2(ifftshift(F_Escat))); % Retrieving scattered beam with inverse Fourier transform.
figure, imagesc(abs(Escat)), title('Escat');

%------- Scattering part due to the cross wave --------%
Ex_response2 = Exy;
F_Ein2=fftshift(fft2(ifftshift(S.*Ein2))); %Fourier transform of the input field.
figure, imagesc(abs(Exy)), title('Exy');
% figure, imagesc(abs(F_Ein));

F_response2 = fftshift(fft2(ifftshift(Ex_response2))); %Fourier transform of the response
% figure, imagesc(abs(F_response));

F_Escat2 = F_Ein2.*F_response2; %Direct scattered field calculation
Escat2=fftshift(ifft2(ifftshift(F_Escat2))); % Retrieving scattered beam with inverse Fourier transform.
% figure, imagesc(abs(Escat));

%-------------------------------------------------------%

%-----------------------Add direct beam------------------%
phz_fft = fft_phase_screen(D, N, lamda1);
phase = exp(1i.*phz_fft);
% phase = 1;
Ein_phase = phase.*Ein;
H2=fft2(fftshift(Ein_phase)); %H2 is the Fourier transform of the input field.
P2=1.*(sqrt(fx.^2+fy.^2)<=n1/lamda1)-1i.*(sqrt(fx.^2+fy.^2)>n1/lamda1);% This step assigns a value of 1 for ?fx2 + fy2 <= n1/lamda1 and of i for ?fx2 + fy2 >n 1/lamda1.
P2=exp(-1i*2*pi*d.*P2.*sqrt(abs((n1/lamda1)^2-fx.^2-fy.^2)));%Transfer function of free space propagation, see Eq. (13).
P2=fftshift(P2); %equation
G2=H2.*P2; % Direct beam calculation executed at the Fourier domain.
Edirect=ifftshift(ifft2(G2)); % Retrieving direct beam with inverse Fourier transform.
% figure, imagesc(abs(Edirect));

%--------------------------------------------------------%

Eout = (Edirect+Escat+Escat2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f1 here uses the Mie scattering library function Mie_S12 to give the required response
% n,cosTh, m, x, u, z, kz are as defined above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Refer to:
%% Lecture 6. Scattering and absorption by aerosol and cloud particles. 
%% Mie theory. Main radiation law (Beer-Bouger- Lambert law) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Efx] = f1(n,cosTh,m,x,u,z,kz)
R = z./cosTh;
% u=1;
S = Mie_S12(m,x,u,1);
S12 = S;
% Efx = exp(-1i*kz*z./cosTh).*S12;
Efx = (exp(-1i*kz*R+1i*kz*z)./(1i*kz*R)).*S12;
end


% function [Efx] = f1(n,cosTh,m,x,u,z,kz)
% R = z./cosTh;
% [m1,n1] = size(u);
% S12 = zeros(size(u));
% for i=1:m1
%     for j=1:n1
%         ud = u(i,j);
%         S = My_Mie_S12(m,x,ud);
%         S12(i,j) = S(n+1);
%     end
% end
% % Efx = exp(-1i*k*z/cosTh)*S2;
% Efx = (exp(-1i*kz*R+1i*kz*z)./(1i*kz*R)).*S12;
% end

