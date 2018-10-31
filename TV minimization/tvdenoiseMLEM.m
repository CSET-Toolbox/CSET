function  u = tvdenoiseMLEM(f, lambda, s, Tol)  
%  function  [u, divp] = tvdenoiseMLEM(f, lambda, s, divp, Tol)  


if nargin < 4
    Tol = 1e-2;
end

N = size(f);
id = [2:N(1),N(1)];
iu = [1,1:N(1)-1];
ir = [2:N(2),N(2)];
il = [1,1:N(2)-1];

p1 = ones(size(f));
p2 = ones(size(f));
divp = ones(size(f));
lastdivp = zeros(size(f));
% auxdivp = zeros(size(f));
% lastauxdivp = zeros(size(f));
% tk = 1;

if length(N) == 2           % TV denoising
   j=0;
   % dt = 0.25/lambda/8;
   aux = s.*f;
   dt = (min(s(:))-4*lambda)^2/(8*lambda*max(aux(:))) ;
   % dt = max(s(:))./(8*lambda*max(f(:)))
   while ((norm(divp(:) - lastdivp(:),inf) > Tol*dt) && (j<2000))
        %z = s.*f./(s+min(max(lambda*divp,-s+10^(-5)),s-10^(-5)));
        z = s.*f./(s+ lambda*divp);
        %z = max(z,f*10^(-5));
        z(z<0) = 0;
        z1 = (z(:,ir) - z);
        z2 = (z(id,:) - z);
        denom = 1+dt*sqrt(z1.^2 + z2.^2);
        p1 = (p1 - dt*z1)./denom;
        p2 = (p2 - dt*z2)./denom;
        
%         lastauxdivp = auxdivp;
%         auxdivp = p1 - p1(:,il) + p2 - p2(iu,:);
        lastdivp = divp;
        divp = p1 - p1(:,il) + p2 - p2(iu,:);
%         tkp1 = (1+sqrt(1+4*tk^2))/2;
%         divp = auxdivp +(tk-1)/tkp1*(auxdivp-lastauxdivp);
        j = j+1;
   end
end

if length(N) == 3           % TV denoising
    ibk = [2:N(3),N(3)];
    ifr = [1,1:N(3)-1];
    p3 = ones(size(f));
   j=0;
   aux = s.*f;
   dt = (min(s(:))-6*lambda)^2/(12*lambda*max(aux(:))) ;
   while ((norm(divp(:) - lastdivp(:),inf) > Tol*dt) && (j<200))
        z = s.*f./(s+ lambda*divp);
        z(z<0) = 0;
        z1 = (z(:,ir,:) - z);
        z2 = (z(id,:,:) - z);
        z3 = (z(:,:,ibk) - z);
        denom = 1 + dt*sqrt(z1.^2 + z2.^2 + z3.^2);
        p1 = (p1 - dt*z1)./denom;
        p2 = (p2 - dt*z2)./denom;
        p3 = (p3 - dt*z3)./denom;
   
        lastdivp = divp;
        divp = p1 - p1(:,il,:) + p2 - p2(iu,:,:) + p3 - p3(:,:,ifr);
        j = j+1;
   end
  
end 

u = max(s.*f./(s + lambda*divp) , 0);
   
end