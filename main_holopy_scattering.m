%% This code can be used as a main file for running the Mie code

% Here,
% L: Length of diffuser
% D: Diffuser width
% deltax: Distance between two sampling points
% N: Number of sampling points
% F: Focal length
% G: Number of particles in LxL space screen
% The unit of distance measurement is micron

%% In this section, you can set all the required values of the parameters

load Python_Data_02_03_19

% L=50; %Length of diffuser
L = 10000; %NEW: 10mm strip
% D = 50; %Diffuser width
D = 2000; %NEW: 2mm far
deltax = 10; %Distance between two sampling points.
N=ceil(D/deltax)+1; % Number of sampling points.
factor = 1;
N1=ceil(D/deltax)*factor+1; % Number of sampling points.
W = D/5;
% F = 100;
F = 10000; %NEW: Set focal length to 10mm
lamda1 = 0.75; % wavelength: 750nm
k = 2*pi/lamda1;
% G=1e6; % Total particle numbers.
G=1000; %NEW: number of particles in LxL space screen
M = 1; %Uncomment this when not using My_Main_Mie_nig
x=linspace(-D/2,D/2,N1); % Spatial position of the diffuser.
[x, y] = meshgrid(x);
map = ones(size(x));
% map = exp(+1i*k.*((x.^2+y.^2)./(2*F))); %Map that acts as if a focusing lens is put in front of the beam

%% Calculation of all the fields

% Create screen with the required number of particles randomly distributed
h1 = My_Step1ScreenParticles(G,N,M,D,deltax,factor);
h = reshape(h1(1,:,:),[N1,N1]);
figure, imagesc(h)
% keyboard;

% Get the field E1 just outside the diffuser
E1 = My_Step12Integrator(1,h1,map,M,D,L,deltax,factor, python_data_x, python_data_y);
figure, imagesc(abs(E1)), title('E1');
% keyboard;

% Get the field Ef1 which is in the far field zone
Ef1 = focusBeam(E1,F,F-L,D,deltax,factor);
figure, imagesc(abs(Ef1)), title('Ef1');
% keyboard;

% Get the field E0 just outside the diffuser
E0 = My_Step12Integrator(0,h1,map,M,D,L,deltax,factor, python_data_x, python_data_y);
figure, imagesc(abs(E0)), title('E0');
% keyboard

% Get the field Ef0 which is in the far field zone
Ef0 = focusBeam(E0,F,F-L,D,deltax,factor);
figure, imagesc(abs(Ef0)), title('Ef0');
% keyboard;

% Get the combined field E just outside the diffuser
E = sqrt(abs(E1).^2+abs(E0).^2);
figure, imagesc(E), title('E');
% keyboard;

% Get the field Ef which is in the far field zone
Ef = sqrt(abs(Ef1).^2+abs(Ef0).^2);
figure, imagesc(abs(Ef)), title('Ef');
% keyboard;