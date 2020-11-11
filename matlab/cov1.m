function[cov_1]=cov1(obx,prex,z)% the convariance k
ob_dim=length(obx);
dim_pre=length(prex);
% par_dim=length(z);% z is an one dimension vactor
cov_1(1:dim_pre,1:ob_dim)=0;
for i=1:dim_pre
    for j=1:ob_dim
        s=(prex(i)-obx(j))^2;
        s=s/(2*z(1)*z(1));
        cov_1(i,j)=z(2)*z(2)*exp(-s);
    end 
end 
end