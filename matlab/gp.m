clc;
clear;
%% initial 
dim_pre=161;%the dimension of prediction
ob_dim=20;%the dimension of calibration
max_gen=80000;
par_dim=3;%the dimension of parameter
%设置不同的超参
% z(1:par_dim)=[1,1,0.1];%parameter
% z(1:par_dim)=[0.3,1.08,0.00005];
z(1:par_dim)=[3.0,1.16,0.89];
in_dim=1;
%% load data
% load('error0.txt');
% load('error1.txt');
% load('error2.txt');
load('ob.txt');%load observation
load('prex.txt');%load the prediction point
%% the covariance
cov_0(:,:)=cov0(ob(:,1),z);%k(x,x)
cov_1(:,:)=cov1(ob(:,1),prex(:),z);%k(x*,x)
cov_2(:,:)=cov2(prex(:),z);%k(x*,x*)
%% ave and variance
av_b=cov_1/(cov_0)*ob(:,2);%average
var_b=cov_2-cov_1/(cov_0)*cov_1';%variance
%% random number
[r,p]=chol(var_b);
var_b=diag(diag(var_b));
b_pre=mvnrnd (av_b,var_b,max_gen);