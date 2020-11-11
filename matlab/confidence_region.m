  clc
  
  %% prediction point
  %load data
%   load('resultoutput.txt');
%   load('b.txt');
  x=prex;
  result=sort(b_pre);%sort 
  low=result(2000,:);%the lower
  upper=result(78000,:);%the upper
  f = [low'; flipdim(upper',1)];
  fill([x; flipdim(x,1)], f, [7 7 7]/8)%[7 7 7]/8É«ºÅ
  hold on; 
% plot(x, mod_ave);
  ave=mean(result);
  plot(x, ave,'r');
  %% calibration point
    %load data
%     load('a.txt');
    plot(ob(:,1), ob(:,2),'+');
    hold on;
%   hold on;
%   %   y1=sin(3*xs);
%   plot(xo(1:35),yo(1:35),'*');
%   hold on;
%   plot(xo(36:end),yo(36:end),'o')