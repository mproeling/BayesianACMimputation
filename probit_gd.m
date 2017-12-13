% PURPOSE: demo of probit_g
%       Gibbs sampling for probit heteroscedastic estimation
%---------------------------------------------------
% USAGE: probit_gd
%---------------------------------------------------

clear all;

rng(2017);

n=101;
k = 3;
evec = randn(n,1);
tt=1:n;

x = randn(n,k);
x(1:n,1) = ones(n,1);

b = ones(k,1);
b(3,1) = -2.0;
b(2,1) = 2.0;
b(1,1) = -0.5;

y = x*b + 0.2*evec;
yc = zeros(n,1);
% now censor the data
for i=1:n
 if y(i,1) > 0
 yc(i,1) = 1;
 else
 yc(i,1) = 0;
 end;
end;

% add outliers
%x(50,2) = 5;
%x(75,3) = 5;

Vnames = strvcat('y','constant','x1','x2');

prior.rval = 40;     % heteroscedastic prior
ndraw = 10100;
nomit = 100;

% original model
result = probit_g(yc,x,ndraw,nomit,prior);
mean(result.bdraw) % -0.6036    2.2556   -2.3779

% model with 1 person imputation
result = probit_g_mp(yc,x,ndraw,nomit,prior);
mean(result.bdraw) % -0.6078    2.2324   -2.3547

h1 = histogram(result.imputation(:,5));
hold on
h2 = histogram(result.imputation(:,6));
hold off
legend show
 
% glmfit(x(:,[2 3]),yc,'binomial', 'link', 'logit')

plot(tt,result.vmean);
title('vi-estimates');
pause;

prt(result,Vnames);

tt=1:n;
[ys yi] = sort(result.y);
plot(tt,ys,tt,result.yhat(yi,1),'--');

plot(tt,result.ymean,tt,y,'o');