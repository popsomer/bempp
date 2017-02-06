% Investigate correctness of the asymptotic expansion of the integral for one obstacle

%% Initialising
clearvars
% close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

obst = 1;
par = getObst(obst);
taus = 0.4; % \tau^*

gamx = zeros(maxOrder,1);
gamy = zeros(maxOrder,1);

if obst == 1
    %'par', @(t) repmat([0;0], 1, length(t)) + [0.5*cos(2*pi*t); 0.5*sin(2*pi*t)]
%     gamx = [cos(2*pi*taus); -2*pi*sin(2*pi*taus); -(2*pi)^2*cos(2*pi*taus); (2*pi)^3*sin(2*pi*taus); (2*pi)^4*cos(2*pi*taus)];
    gamx = [cos(2*pi*taus); -2*pi*sin(2*pi*taus); -(2*pi)^2*cos(2*pi*taus)/2; (2*pi)^3*sin(2*pi*taus)/2/3; (2*pi)^4*cos(2*pi*taus)/2/3/4];
    gamx = gamx*0.5;
    gamy = [sin(2*pi*taus); 2*pi*cos(2*pi*taus); -(2*pi)^2*sin(2*pi*taus)/2; -(2*pi)^3*cos(2*pi*taus)/6; (2*pi)^4*sin(2*pi*taus)/24];
    gamy = gamy*0.5;
end
binom = @(x,n) prod(x:-1:(x-n+1) )/prod(1:n)*0^((n<0) + abs(rem(n,1)) ); % Not to be used vectorized (and n should be a positive integer)
% maxOrder = 4;
maxOrder = length(gamx)-1;

%% Check series expansion of Jacobian
% Omx = nan*gamx;
Omx = nan*gamx(1:end-1);
Omy = nan*gamy(1:end-1);
% zeta = zeros(maxOrder,maxOrder); % column index corresponds to power of y
zeta = nan*zeros(maxOrder-1,maxOrder-1); % column index corresponds to power of y
%zeta(1,1) = 1;
for j = 1:maxOrder
    Omx(j) = sum(transpose((1:j).*(j:-1:1)).*gamx(2:(j+1)).*gamx((j+1):-1:2));
    Omy(j) = sum(transpose((1:j).*(j:-1:1)).*gamy(2:(j+1)).*gamy((j+1):-1:2));
%     Omx(j) = sum((1:j).*gamy(2:(j+1)).*(j:-1:1).*gamy((j+1):-1:2));
end
for j=1:maxOrder-1
    zeta(1,j) = (Omx(j+1)+Omy(j+1))/(Omx(1)+Omy(1));
%     zeta(j,1) = (Omx(j+1)+Omy(j+1))/(Omx(1)+Omy(1));
end

for n=2:maxOrder
    for i=1:maxOrder-1
        zeta(n,i) = sum(zeta(1,1:i).*zeta(n-1,i:-1:1) );
%         zeta(i,n) = sum(zeta(1:i,1).*zeta(i:-1:1,n-1) );
    end
end
% d = ones(size(gamx));
d = zeros(size(gamx));
d(1) = 1;
for n = 2:maxOrder
    for i=1:(n-1)%n
%         d(n) = d(n) + binom(1/2,n-i)*zeta(i,n-i);
        d(n) = d(n) + binom(1/2,n-i)*zeta(n-i,i);
    end
end
d = d*sqrt(Omx(1)+Omy(1)) %Should be [pi; 0; 0; ...] for obst=1

%% Series of the distance function
Lambdax = nan*gamx(1:end-1);
Lambday = nan*gamx(1:end-1);
for j=1:maxOrder
    Lambdax(j) = sum(gamx(2:(j+1)).*gamx((j+1):-1:2));
    Lambday(j) = sum(gamy(2:(j+1)).*gamy((j+1):-1:2));
end
% \Lambda_j^{x/y} & = \sum_{i=1}^j \Gamma_{i+1}^{x/y} \Gamma_{j-i+2}^{x/y} \\
z = nan*zeros(maxOrder-1,maxOrder-1); % column index corresponds to power of y
for i=1:maxOrder-1
    z(i,1) = (Lambdax(i+1)+Lambday(i+1))/(Lambdax(1)+Lambday(1));
end
% z_i^1 & = \frac{\Lambda_{i+1}^x + \Lambda_{i+1}^y}{\Lambda_1^x + \Lambda_1^y} \\
for n=2:maxOrder
    for i=1:maxOrder-1
        z(i,n) = sum(z(1:i,1).*z(i:-1:1,n-1) );
    end
end
% z_i^n & = \sum_{k=1}^i z_k^1 z_{i-k+1}^{n-1} \\
f = zeros(size(gamx));
f(1) = 1;
for n = 2:maxOrder
    for i=1:(n-1)%n
        f(n) = f(n) + binom(1/2,n-i)*zeta(i,n-i);
    end
end
f = f*sqrt(Lambdax(1) + Lambday(1))
% f_n & = \sqrt{\Lambda_1^x + \Lambda_1^y}\sum_{i=1}^{n-1} {1/2 \choose n-i} z_i^{n-i} \quad n \geq 2
