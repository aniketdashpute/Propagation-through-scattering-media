%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code creates M screens each with G/M particles 
% with random distribution of particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h2 = My_Step1ScreenParticles(G,N,M,D,deltax,factor)

% keyboard;
N1 = ceil(D/deltax)*factor+1; % Number of sampling points.
g=ceil(G/M); % Particle numbers per screen.
% The procedures below generate a 2D random position map with g particles.
h1=zeros(M,N,N); % CHANGE
for i1=1:M
    for m=1:g
        b=ceil(N*rand(1,1));
        c=ceil(N*rand(1,1));
        h1(i1,b,c)=1;
    end
end
h2=h1;
% %-------------------------------%
% keyboard;
% h=reshape(h1,[N,N]);
% figure, imagesc(h);
% h2 = zeros(M,N1,N1);
% for i1=1:M
%     for i=1:N1
%         for j=1:N1
%             a=ceil(i/factor);
%             b=ceil(j/factor);
%             h2(i1,i,j) = h1(1,a,b);
%         end
%     end
% end
% %-------------------------------%
% h=reshape(h2,[N1,N1]);
% figure, imagesc(h);
% keyboard;
%%--------------------------------------------- %%