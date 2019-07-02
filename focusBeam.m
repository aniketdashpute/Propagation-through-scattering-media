function out =  focusBeam(Ein,F,d,D,deltax,factor)
% number = 0;
n1=1; % Refractive index of the medium.
lamda1 = 0.5;

N=ceil(D/deltax)*factor+1; % Number of sampling points.

E = Ein;

fx=1/deltax*(-0.5:1/(N-1):0.5)*factor; % Spatial frequencies.
[fx, fy] = meshgrid(fx);

H2=fft2(fftshift(E)); %H2 is the Fourier transform of the input field.

P2=1.*(sqrt(fx.^2+fy.^2)<=n1/lamda1)-1i.*(sqrt(fx.^2+fy.^2)>n1/lamda1);% This step assigns a value of 1 for ?fx2 + fy2 <= n1/lamda1 and of i for ?fx2 + fy2 >n 1/lamda1.
P2=exp(-1i*2*pi*d.*P2.*sqrt(abs((n1/lamda1)^2-fx.^2-fy.^2)));%Transfer function of free space propagation, see Eq. (13).
P2=fftshift(P2); %equation.

G2=H2.*P2; % Direct beam calculation executed at the Fourier domain.
g2=ifftshift(ifft2(G2)); % Retrieving direct beam with inverse Fourier transform.

out = g2;
% figure, imshow(abs(out)),title(d);
end