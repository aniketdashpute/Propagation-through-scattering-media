%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes care of returning phase screens given D, N, lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phz_fft = fft_phase_screen(D, N, lambda)

% D = D*(10^-3);
% lambda = lambda*(10^-3);

del_x = D/N;
[x,y]= meshgrid(-D/2:del_x:D/2-del_x);
B=1/del_x;
del_f = 1/(N*del_x);                      % frequency grid spacing [1/m]
[fx,fy]=meshgrid(-B/2:del_f:B/2-del_f);   % frequency grid [1/m]
[th, f]=cart2pol(fx,fy);
% keyboard;

% fm = 5.92/lo/(2*pi);                      % inner scale frequency [1/m]
% fo=1/Lo;                                  % outer scale frequency [1/m]

% constant = 0.001;
% refidx=0.033.*Cn2.*(f.^2+fo.^2).^(-11/6).*(2*pi).^(-11/3);   % modified Von Karman
%refidx=0.033.*Cn2.*(f.^2+fo.^2).^(-11/6).*(2*pi).^(-11/3).*exp(-(f*lo).^2);   % Von Karman
% refidx=0.033.*10^(-15).*(f).^(-2.5);           % Kolmogorovxxx

%% ---------------Kolmogorov-------------------- %%
k=(2*pi)/lambda;                          % wave number
del_z = 5e-3;
constant=2.*pi.*k.^2.*del_z;
Cn2 = 10^(-4);
refidx=0.033.*Cn2.*(f).^(-11/3).*(2*pi).^(-11/3);           % Kolmogorov
PSD_phi=constant.*refidx;
PSD_phi(ceil(N/2+1),ceil(N/2+1)) = 0;      

%% -------------------SVR----------------------- %%

% f = f.*(10^3);
% 
% ff = reshape(f,[N*N,1]);
% ff = log(ff);
% csvwrite('ff.csv',ff);
% keyboard;

% 
% load SVR_Results;
% f_predict = reshape(ff_predict,[N,N]);
% PSD_phi=f_predict;
% PSD_phi(ceil(N/2+1),ceil(N/2+1)) = 0;                                   %putting middle pixel equal to zero

%% ---------------------------------------------- %%
A = 2*pi*del_f.*(randn(N) + 1i*randn(N)).* sqrt(PSD_phi);   %Filtering with required spectra
phz1 = fftshift(ifft2(ifftshift(A))).*N.*N;
phz_fft=real(phz1);

end