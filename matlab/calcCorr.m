% Calculate the correlations.
% Input
%   par       - The structure containing k,par, and so on
%   c         - The solution vector to compute the correlations with
%   Tcor      - The width of the windows
%   pd        - The percentage of support of the windows with a decay
%   [tim      - If present and nonempty, tim(1) is the number of seconds after which 
%                  to print tim(2) + expected days for this method, assuming tic has been started]
%   [A        - If present, the system matrix to compute the correlations with]
% Output
%   R         - The correlations where each row corresponds to par.colltau
%   sigm      - The locations of the windows, corresponding the columns of R
%   [obbounds - Additional information for multiple scattering obstacles]
function [R, sigma, obbounds] = calcCorr(par, c, Tcor, pd, tim, A)
if ~exist('tim','var') || isempty(tim)
    tim = [inf, 0];
end
fr = 1.5; % factor of number of columns over rows of R
if isfield(par, 'obsts')
    sigma = [];
    obbounds = zeros(2,length(par.obsts));
    for obst = 1:length(par.obsts)
        obbounds(1,obst) = length(sigma)+1;
        sigma = [sigma linspace(0,1,round(par.obsts(obst).N*fr) )];
        obbounds(2,obst) = length(sigma);
    end
    sr = length(sigma);
    R = zeros(par.N,sr);
else
    R = zeros(par.N,round(par.N*fr));
    sigma = linspace(0,1,size(R,2) );
end
sr = size(R,2);
prevToc = toc;

if exist('A','var') && isfield(par, 'obsts')
    % We go over the columns of R so that we can reuse j1, j2 and wi
    for roi = 1:sr
        if (toc-prevToc > tim(1) )
            prevToc = toc;            
            display(['k=' num2str(par.k) ': ' num2str(roi/sr,'%7.3f')  '=R%, now=' datestr(now) ', est. # sec. left for R=' ...
                num2str(toc*(sr-roi)/roi) ', est. end ' datestr(tim(2)+prevToc*sr/roi/24/3600)])
        end
        obsin = find((roi >= obbounds(1,:)) & (roi <= obbounds(2,:)));
        T1 = sigma(roi)-Tcor;
        T2 = sigma(roi)+Tcor;
        [j1loc, j2loc, noSplit] = bounds2Ind(par.obsts(obsin).colltau,T1,T2);
        j1 = j1loc+par.r(1,obsin)-1;
        j2 = j2loc+par.r(1,obsin)-1;
        if noSplit
            wi = chi(par.obsts(obsin).colltau(j1loc:j2loc),T1, T1+pd*Tcor, T2-pd*Tcor, T2,0);
            for i = 1:par.N
                R(i,roi) = (wi.*A(i,j1:j2))*c(j1:j2);
            end
        else
            wi1 = chi(par.obsts(obsin).colltau(j1loc:par.obsts(obsin).N),T1, T1+pd*Tcor, T2-pd*Tcor, T2,1);
            wi2 = chi(par.obsts(obsin).colltau(1:j2loc),T1, T1+pd*Tcor, T2-pd*Tcor, T2,1);
            for i = 1:par.N
                R(i,roi) = (wi1.*A(i,j1:par.r(2,obsin)))*c(j1:par.r(2,obsin)) + (wi2.*A(i,par.r(1,obsin):j2))*c(par.r(1,obsin):j2);
            end
        end
    end
    return
elseif exist('A','var')
    % We go over the columns of R so that we can reuse j1, j2 and wi
    for roi = 1:sr
        if (toc-prevToc > tim(1) )
            prevToc = toc;
            display(['k=' num2str(par.k) ': ' num2str(roi/sr,'%7.3f')  '=R%, now=' datestr(now) ', est. # sec. left for R=' ...
                num2str(toc*(sr-roi)/roi) ', est. end ' datestr(tim(2)+prevToc*sr/roi/24/3600)])
        end
        T1 = sigma(roi)-Tcor;
        T2 = sigma(roi)+Tcor;
        [j1, j2, noSplit] = bounds2Ind(par.colltau,T1,T2);
        if noSplit
            wi = chi(par.colltau(j1:j2),T1, T1+pd*Tcor, T2-pd*Tcor, T2,0);
            for i = 1:par.N
                R(i,roi) = (wi.*A(i,j1:j2))*c(j1:j2);
            end
        else
            wi1 = chi(par.colltau(j1:par.N),T1, T1+pd*Tcor, T2-pd*Tcor, T2,1);
            wi2 = chi(par.colltau(1:j2),T1, T1+pd*Tcor, T2-pd*Tcor, T2,1);
            for i = 1:par.N
                R(i,roi) = (wi1.*A(i,j1:par.N))*c(j1:par.N) + (wi2.*A(i,1:j2))*c(1:j2);
            end
        end
    end
    return
end

% Compute the correlations as an integral through a simple Riemann-sum without normalisation.
% We could compute R column by column because then the window around sigma(roi) would be constant, so that we do not have to
% recompute the same window indices j1 and j2 for each row. But this would need to recompute the integrand with the
% (expensive) Bessel function for each column.
if ~isfield(par, 'obsts') % Single scattering obstacle without using A
    u = [par.t(end-par.dbf:end-1)-1, par.t, par.t(2:par.dbf+1)+1]; % the extended knots
    c1ip = deboor(par.dbf, u, [c; c(1:par.dbf)], sigma);
    
    % Find the intervals of u where the respective sigma lies.
    for i = 1:par.N
        if (toc-prevToc > tim(1) )
            prevToc = toc;
            display(['k=' num2str(par.k) ': ' num2str(i/par.N,'%7.3f')  '=R%, now=' datestr(now) ', est. # sec. left for R=' ...
                num2str(toc*(par.N-i)/i) ', est. end ' datestr(tim(2)+prevToc*par.N/i/24/3600)])
        end
        tc = par.colltau(i);
        integrand = 1i/4.*besselh(0, 1, par.k*sqrt(sum((par.par(tc*ones(size(sigma))) ...
            -par.par(sigma) ).^2, 1)) ).*c1ip.*par.gradnorm(sigma);
        integrand(isnan(integrand)) = 0;
        for roi = 1:size(R,2)
            [j1, j2, noSplit] = bounds2Ind(sigma,sigma(roi)-Tcor,sigma(roi)+Tcor);
            if noSplit
                wi = chi(sigma(j1:j2),sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,0);
                R(i,roi) = sum(wi.*integrand(j1:j2) );
            else
                wi1 = chi(sigma(j1:sr),sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,1);
                wi2 = chi(sigma(1:j2), sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,1);
                R(i,roi) = sum(wi1.*integrand(j1:sr) ) + sum(wi2.*integrand(1:j2) );
            end
        end
    end
end
% Multiple scattering obstacle without using A
c1ip = NaN*sigma;
for obst = 1:length(par.obsts)
    u = [par.obsts(obst).t(end-par.obsts(obst).dbf:end-1)-1, par.obsts(obst).t, par.obsts(obst).t(2:par.obsts(obst).dbf+1)+1]; % the extended knots
    c1ip(obbounds(1,obst):obbounds(2,obst)) = deboor(par.obsts(obst).dbf, u, ...
        [c(par.r(1,obst):par.r(2,obst)); c(par.r(1,obst):(par.r(1,obst) +par.obsts(obst).dbf-1))], sigma(obbounds(1,obst):obbounds(2,obst)) );
end
for i = 1:par.N
    if (toc-prevToc > tim(1) )
        prevToc = toc;
            display(['k=' num2str(par.k) ': ' num2str(i/par.N,'%7.3f')  '=R%, now=' datestr(now) ', est. # sec. left for R=' ...
                num2str(toc*(par.N-i)/i) ', est. end ' datestr(tim(2)+prevToc*par.N/i/24/3600)])
    end
    obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
    tc = par.obsts(obsin).colltau(i-par.r(1,obsin)+1);
    colx = par.obsts(obsin).par(tc);
    for obst = 1:length(par.obsts)
        cusigma = sigma(obbounds(1,obst):obbounds(2,obst));
        integrand = 1i/4.*besselh(0, 1, par.k*sqrt(sum((repmat(colx,1,length(cusigma)) -par.obsts(obst).par(cusigma) ).^2, 1)) ).*...
            c1ip(obbounds(1,obst):obbounds(2,obst)).*par.obsts(obsin).gradnorm(cusigma);
        integrand(isnan(integrand)) = 0;
        for roi = obbounds(1,obst):obbounds(2,obst)
            [j1, j2, noSplit] = bounds2Ind(cusigma,sigma(roi)-Tcor,sigma(roi)+Tcor);
            if noSplit
                wi = chi(cusigma(j1:j2),sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,0);
                R(i,roi) = sum(wi.*integrand(j1:j2) );
            else
                wi1 = chi(cusigma(j1:length(cusigma)),sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,1);
                wi2 = chi(cusigma(1:j2), sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,1);
                R(i,roi) = sum(wi1.*integrand(j1:length(cusigma)) ) + sum(wi2.*integrand(1:j2) );
            end
        end
    end
end

end

% Evaluate the spline of degree k determined by the knots t(n+2*k+1) and coefficients c(n+k) at the points in x.
function s = deboor(k, t, c, x)

s = zeros(size(x));
for in=1:length(x)
    % Determine the interval where x(in) lies
    if x(in) == t(length(t)-k)
        j = length(t) - 2*k - 2;
    else
        j = find(x(in) >= t, 1, 'last') - k - 1;
    end
    d=c(j-k+k+1:j+k+1);
    for r=1:k
        for i=j:-1:j-k+r
            idx = i+k+1;
            alfa = (x(in) - t(idx)) / (t(idx+k+1-r)-t(idx));
            d(idx-j) = alfa*d(idx-j) + (1-alfa)*d(idx-j-1);
        end
    end
    s(in) = d(k+1);
end

end
