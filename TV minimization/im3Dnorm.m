function n=im3Dnorm(img,normind,varargin)
% IMAGE3DNORM computes the desired image norm
%   IMAGE3DNORM(IMG,NORMIND) computes the norm if image IMG using the norm
%   defined in NORMING
%
%   IMG         A 3D image
%   NORMIND     'LN': N norm
%               'TV': TV norm
% 


if strcmp(normind(1),'L')
    Nnorm=str2double(normind(2:end));
    n=norm(img(:),Nnorm);
    return;
end

if strcmp(normind,'TV')
    typeGrad=varargin{1};
    if strcmpi(typeGrad,'central')
        [gx,gy,gz]=gradient(img);
    end
    if strcmpi(typeGrad,'forward')
        gx=diff(img,1,1);
        gy=diff(img,1,2);
        gz=diff(img,1,3);
        gx=cat(1,gx,zeros(size(gx(end,:,:))));
        gy=cat(2,gy,zeros(size(gy(:,end,:))));
        gz=cat(3,gz,zeros(size(gz(:,:,end))));

    end
    if strcmpi(typeGrad,'backward')
        gx=diff(img,1,1);
        gy=diff(img,1,2);
        gz=diff(img,1,3);
        gx=cat(1,zeros(size(gx(1,:,:))),gx);
        gy=cat(2,zeros(size(gy(:,1,:))),gy);
        gz=cat(3,zeros(size(gz(:,:,1))),gz);
    end

    g=sqrt(gx.^2+gy.^2+gz.^2); clear gx gy gz;
    n=sum(g(:));
    return;
end


end