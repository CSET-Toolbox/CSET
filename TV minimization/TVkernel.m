%Function: TV minimization 
function In = TVkernel(I, alpha, d, Ntv)
eps = 1e-8;     % in case of zero division
for i = 1:Ntv
    Ipad = padarray(I,[1 1 1]);
    v1 = ((I - Ipad(1:end-2,2:end-1,2:end-1)) + (I - Ipad(2:end-1,1:end-2,2:end-1)) + (I - Ipad(2:end-1,2:end-1,1:end-2)))... 
        ./sqrt((I - Ipad(1:end-2,2:end-1,2:end-1)).^2 + (I - Ipad(2:end-1,1:end-2,2:end-1)).^2 + (I - Ipad(2:end-1,2:end-1,1:end-2)).^2 + eps);
    v2 = (Ipad(3:end,2:end-1,2:end-1)-I)./sqrt((Ipad(3:end,2:end-1,2:end-1) - I).^2 + (Ipad(3:end,2:end-1,2:end-1) - Ipad(3:end,1:end-2,2:end-1)).^2 + (Ipad(3:end,2:end-1,2:end-1) - Ipad(3:end,2:end-1,1:end-2)).^2 + eps);
    v3 = (Ipad(2:end-1,3:end,2:end-1)-I)./sqrt((Ipad(2:end-1,3:end,2:end-1) - Ipad(1:end-2,3:end,2:end-1)).^2 + (Ipad(2:end-1,3:end,2:end-1) - I).^2 + (Ipad(2:end-1,3:end,2:end-1) - Ipad(2:end-1,3:end,1:end-2)).^2 + eps);
    v4 = (Ipad(2:end-1,2:end-1,3:end)-I)./sqrt((Ipad(2:end-1,2:end-1,3:end) - Ipad(1:end-2,2:end-1,3:end)).^2 + (Ipad(2:end-1,2:end-1,3:end) - Ipad(2:end-1,1:end-2,3:end)).^2 + (Ipad(2:end-1,2:end-1,3:end) - I).^2 + eps);
    v = v1 - v2 - v3 -v4;
    
    norm3d = im3Dnorm(v,'L2');
    
    v = v / norm3d;
    I = I - alpha*d*v;
end

In = I;

end

