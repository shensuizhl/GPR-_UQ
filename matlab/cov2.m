function[cov_2]=cov2(prex,z)% the convariance k
dim_pre=length(prex);
% par_dim=length(z);% z is an one dimension vactor
cov_2(1:dim_pre,1:dim_pre)=0;
for i=1:dim_pre
    for j=1:dim_pre
        s=(prex(i)-prex(j))^2;
        s=s/(2*z(1)*z(1));
        cov_2(i,j)=z(2)*z(2)*exp(-s);
    end 
end 
end