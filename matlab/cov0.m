function[cov_0]=cov0(obx,z)
ob_dim=length(obx);
% par_dim=length(z);% z is an one dimension vactor
cov_0(1:ob_dim,1:ob_dim)=0;
for i=1:ob_dim
    for j=1:ob_dim
        s=(obx(i)-obx(j))^2;
        s=s/(2*z(1)*z(1));
        if (i==j)
            cov_0(i,j)=z(2)*z(2)*exp(-s)+z(3)*z(3);
        else 
            cov_0(i,j)=z(2)*z(2)*exp(-s);
        end
    end 
end 
end