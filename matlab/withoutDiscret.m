% Use a direct computation of the integrals using the exact solution without discretizing the boundary to show one can use window functions.

function withoutDiscret() % This is actually a script but we need short functions in it.
%% Initialising
percDecay = 0.5; % Percentage of the window for the C-inf decay: 0 means a block window and 1 means not identically one on any interval.

par = getObst(0);
x0 = par.par(0);
a = x0(1); % Radius
tol = 1e-8; % Tolerance for computing quadgk

ks = 2.^(6:15);
Ts = linspace(0.04,0.4,25);
shad = 0.91; % Parametrization of the point in the shadow region
pts = [par.par(0.56), par.par(0.25), par.par(shad), [-1.1; -0.4]];

%% Exact solutions for circle, from ietest/highf/circle_field_solution.m
% Input
%	pt - The (possibly internal/external) point
% Output
%	q  - The scattered field
function q = scatField(pt)
        r = sqrt(pt(1)^2+pt(2)^2);
        theta = atan2(pt(2), pt(1));
        if theta < 0
            theta = theta + 2*pi;
        end
        ost = ones(size(theta));
        q2 = 0*ost;
        q = besselj(0,par.k*a*ost)./besselh(0,1,par.k*a*ost).*besselh(0,1,par.k*r);
        % We could check norm of the term added but this might be zero for odd n at theta = pi/2
        for n = 1:1.7*par.k+50
            q1 = q;
            q = q + 2*(1i)^n*besselj(n,par.k*a*ost)./besselh(n,1,par.k*a*ost).*besselh(n,1,par.k*r).*cos(n*theta);
            if norm(q-q2) < tol
                break; % No significant changes in last two iterations so quit
            end
            q2 = q1;
        end
        if (norm(2*(1i)^n*besselj(n,par.k*a*ost)./besselh(n,1,par.k*a*ost).*besselh(n,1,par.k*r).*cos(n*theta) ) > tol) || isnan(norm(q) )
            warning([num2str(n) ' = n, Accuracy for scattered field not reached at k = ' num2str(par.k) ]);
        end
end

% Input
%	tau   - Parameter of the boundary
%	norPt - Norm of the (possibly internal/external) point
% Output
%	q     - The density, which is the normal derivative of the scattered (minus incident) field
function q = dens(tau,norPt)
        ost = ones(size(tau));
        q = besselj(0,par.k*a*ost)./besselh(0,1,par.k*a*ost).*par.k.*besselh(-1,1,par.k*a*ost) ;
        if norPt > a*(1-10*eps)
            % The evaluation point is not inside the circle, so add the normal derivative of the incident plane wave.
            q = q -1i*par.k*cos(2*pi*tau).*exp(1i*par.k*a*cos(2*pi*tau));
        end
        q2 = 0*ost;
        for n = 1:1.7*par.k+50
            q1 = q;
            q = q + 2*(1i)^n*besselj(n,par.k*a*ost)./besselh(n,1,par.k*a*ost).*par.k.*...
                (besselh(n-1,1,par.k*a*ost) -n/(par.k*a)*besselh(n,1,par.k*a*ost)).*cos(2*pi*n*tau);
            if norm(q-q2) < tol
                break; % No significant changes in last two iterations so quit
            end
            q2 = q1;
        end
        if (norm(2*(1i)^n*besselj(n,par.k*a*ost)./besselh(n,1,par.k*a*ost).*par.k.*...
                (besselh(n-1,1,par.k*a*ost) -n/(par.k*a)*besselh(n,1,par.k*a*ost) ).*cos(2*pi*n*tau) ) > tol) || isnan(norm(q) )
            warning([num2str(n) ' = n, Accuracy for density not reached at k = ' num2str(par.k) ]);
        end
end

%% Computing integrals
kl = length(ks);
Tl = length(Ts);
if(max(Ts) > 0.5)
    warning('Not corrected for windows overlapping periodically with themselves');
end
pl = size(pts,2);
% Find the taus giving point on circle closest to pts(pti), but switch to its SP if that point lies in the shadow.
cltau = NaN*ones(pl,1);
for pti = 1:pl
    cltau(pti) = fminbnd(@(tau) sqrt(sum( (repmat(pts(:,pti),1,length(tau))-par.par(tau)).^2,1)),0,1);
end
cltau(3) = 1.5-shad; % SP for collocation point in the shadow

errs = NaN*ones(kl,Tl,pl);
exVal = NaN*ones(kl,pl);
start = now;
for ki = 1:kl
    par.k = ks(ki);
    for pti = 1:pl
        display([datestr(now) '=now, ' num2str(ks(ki)) ' =k, pti = ' num2str(pti)]);
        point = pts(:,pti);
        
        exVal(ki,pti) = - scatField(point); % Take minus so that -BC-sum = 0 and integral-(-sum)=0
        if (pti <= 3) && ( abs((exVal(ki,pti) - par.bc(ks(ki),point))/par.bc(ks(ki),point)) > 10*tol) % exScatField not valid for interior
            warning([num2str(pti) ',' num2str(exVal(ki,pti)) ' = exact solution on boundary is not equal to BC = ' num2str(par.bc(ks(ki),point))])
        end
        integrEx = quadgk(@(tau) 1i/4*besselh(0,1,par.k*sqrt(sum( (repmat(point,1,length(tau))-par.par(tau)).^2,1)))...
            .*dens(tau,norm(point)).*par.gradnorm(tau),0,1, 'RelTol', tol);
        if abs( (integrEx-exVal(ki,pti))/exVal(ki,pti) ) > 10*tol % integral=-sum, also in exterior
            warning([num2str(pti) ',' num2str(exVal(ki,pti)) ' = exact solution in field is not equal to integral = ' num2str(integrEx)])
        end
        for ti = 1:Tl
            errs(ki,ti,pti) = (quadgk(@(tau) 1i/4*besselh(0,1,par.k*sqrt(sum( (repmat(point,1,length(tau))-par.par(tau)).^2,1))) ...
                .*chi(tau,cltau(pti)-Ts(ti), cltau(pti)-Ts(ti)*percDecay, cltau(pti)+ Ts(ti)*percDecay, cltau(pti) + Ts(ti), 1) ...
                .*dens(tau,norm(point)).*par.gradnorm(tau),0,1, 'RelTol', tol)  - exVal(ki,pti))/exVal(ki,pti);
        end
        save withoutDiscret.mat % Save results when running this long computation
    end
    display([datestr(now) '=now, k=' num2str(par.k) ', est. end=' datestr(start + (now-start)*sum(ks.^2)/sum(ks(1:ki).^2) )]);
end


%% Plot the results and find least squares fit T = c k^e to get an approximately fixed error
load withoutDiscret.mat
myTol = 1e-5; % The error for the least squares fit
klsq = zeros(kl,pl);
Tlsq = zeros(kl,pl);
kiLsq = zeros(kl,pl);
locLsq = zeros(kl,pl);
for pti = 1:pl
    for ti = 1:Tl % set err(ki:,ti,pti) to NaN or zero if increases
        [~, loc] = min(errs(:,ti,pti));
        if loc < kl
            errs((loc+1):end,ti,pti) = NaN;
        end
    end
    for ki = 1:kl
        loc = find(abs(errs(ki,:,pti)) < myTol, 1, 'first');
        if (~isempty(loc)) && abs(errs(ki,loc,pti)) > myTol/10
            klsq(ki,pti) = ks(ki);
            Tlsq(ki,pti) = Ts(loc);
            kiLsq(ki,pti) = ki;
            locLsq(ki,pti) = loc;
        end
    end
end

legends = cell(kl,1);
for i = 1:kl
    legends{i} = ['k = ' num2str(ks(i))];
end
set(0,'DefaultFigureWindowStyle', 'docked');
titl = {'Illuminated region', 'Transition region', 'Shadow region', 'Exterior region'};
s = NaN*ones(2,pl); % least squares parameters
linest = {'-','-.','--',':'};
cols = 'rgbmkcy';
fs = 20;
for pti = 1:pl
    figure; set(axes,'LineStyleOrder',linest'); 
    for ki = 1:kl
        semilogy(Ts,abs(errs(ki,:,pti)), [cols(mod(ki,length(cols))+1) linest{mod(ki+3,length(linest))+1}], 'LineWidth', 4);
        hold on;
    end 
    nzk = nonzeros(klsq(:,pti));
    adenkl = abs(diag(errs(nonzeros(kiLsq(:,pti) ), nonzeros(locLsq(:,pti))',pti)));
    semilogy(nonzeros(Tlsq(:,pti)), adenkl, 'k*', 'MarkerSize', 14); 
    xlabel('T'); ylabel('Relative error'); 
    legend(legends(1:10-4*(pti==1)), 'FontSize',fs);
    s(:,pti) = [ones(length(nzk),1), log(nzk)] \ log(nonzeros(Tlsq(:,pti)));
    set(gca, 'FontSize', fs);
end
% Now make the least squares plot
figure;
fs = 28;
hold on; 
% Using ax = gca; ax.ColorOrderIndex = 1; does not make the next plot use the same colors, so Matlab forces us to write loops with 
% explicit colors if we want to combine them with markers or linestyles.
marks = {'v', '+', 'o', 'x', '*', 'h', 'd'};
for pti = 1:pl
     loglog(ks, exp(repmat(s(1,pti),kl,1) + repmat(s(2,pti),kl,1).*log(ks')), [cols(pti) linest{pti}], 'LineWidth', 5);
end
for pti = 1:pl
     loglog(nonzeros(klsq(:,pti)),nonzeros(Tlsq(:,pti)),[marks{pti} cols(pti)], 'MarkerSize', 12);
end
xlabel('k', 'FontSize', fs); legend(titl, 'FontSize', fs); ylabel('T(k)', 'FontSize', fs);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', fs);
s % Determines the least squares functions

end % function withoutDiscret
