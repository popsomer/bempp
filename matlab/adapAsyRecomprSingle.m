% Adaptive Asymptotic Recompression for single scattering obstacles with general geometries: 
% Recompute all l, lambda, rho and r at each wavenumber in an automatic and very efficient way.

%% Initialisation
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

% for R_{m,n}: zeta(tau-sigma) = chi(tau-sigma, -cfT/k, (percDecay-1)*cft/k, (1-percDecay)*cfT/k, cft/k)
percDecay = 1; cfT = 0.8;
% For \tilde{A} around Green singularities: w(t_i,tau) = chi(tau-t_i, -cfTg/k, (pDg-1)*cftg/k, (1-pDg)*cfTg/k, cftg/k)
% pDg = 0.8; cfTg = 1.2;
pDg = 0.5; cfTg = 1.2;
% Added decaying parts on both sides are min(Twi, max(pDg*cfTg/par.k,pbd times the width of the constant part of the window))
pbd = 0.2; Twi = 0.02;

% ksR = [2^4, 2^6]; % The wavenumbers at which to recompute R using A2
% ksR = [2.^(4:6), 104]; % The wavenumbers at which to recompute R using A2
% ksR = [2.^(4:10), (2^11-2^7)]; % The wavenumbers at which to recompute R using A2
% ksR = [2.^(4:7), (2^8-2^5)];
% ksR = 2^9;
ksR = 2.^(4:7);
dk = 8; % Distance between consecutive wavenumbers
% obsts = [1 2 3 4 8];
obsts = 8; %3;

% ks is a list of all wavenumbers and must start with ksR(1:2), then the wavenumbers where we apply
% interpolation between those, then ksR(3) for interpolation of R at ksR(2) and ks(3) and so on
% ks = [ksR, (ksR(1)+dk):dk:(ksR(2)-dk)];
ks = ksR(1);
for a = 2:length(ksR)
    ks = [ks, ksR(a), (ksR(a-1)+dk):dk:(ksR(a)-dk)];
end
% ks = [ks, (ksR(end)+dk):dk:(2*ksR(end))]; % Assuming 2*ksR(end) is the final wavenr
% ks = [ks, (ksR(end)+dk):dk:80] % Assuming this is the final wavenr and steps of dk reach it
% ks = [ks, (ksR(end)+dk):dk:2^11] % Assuming 128 is the final wavenr
% ks = [ksR(1:2), 24, ksR(3), 48, 128];
% ks = sort([ksR(1:2), 24, ksR(3), 48, 128]);
ks = [ksR(1), 24, ksR(2), 48, ksR(3), 88, 104, ksR(4), 160];
kl = length(ks);
maxob = length(obsts);
mti = 2; % Number of tests using iterative solvers, can be zero
avm = 100; % Number of random taus to average the residue of the integral equation over
v = struct('conds', zeros(maxob*kl,2), 'mti', mti, 'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(maxob*kl,2+mti),...
    'perc', zeros(maxob*kl,2), 'errSol', zeros(maxob*kl,2+mti), 'compresErr', zeros(maxob*kl,2),  'errInt', zeros(maxob*kl,2+mti), ...
    'timeSol', zeros(maxob*kl,2+mti), 'nbIter', zeros(maxob*kl,mti), 'timeA', zeros(maxob*kl,4), 'ks', ks);

printtoc = 0.4; % Print info after this number of seconds
powTime = sqrt(2); %1.5; % The estimated complexity for printing the estimated remaining time
ppw = zeros(maxob,1); % Number of points per wavelength to estimate the remaining time
for oi = 1:length(obsts)
    par = getObst(obsts(oi));
    ppw(oi) = par.ppw;
    if (ppw(oi)*min(percDecay*cfT, pDg*cfTg) < 3)
        error(['Decay of the window around Green singularity or sigma_n will have less than three points for obstacle ' num2str(obsts(oi))]);
    end
end
% Following warning not nec if 0.02*16/ks(ki) iso pDg*...
% if Twi < pDg*cfTg/min(ks)
%     warning('Increase Twi')
% end


%% Computation
startAll = now;
expectedEnd = startAll;
% figure; % For plotting max of corr per row for each ksR
for oi = 1:length(obsts)
    obstacle = obsts(oi);
    par = getObst(obstacle);
%     par.xi = par.xi/2; % Test to increase accuracy
    
    for ki = 1:kl
        startk = now;
        idx = (oi-1)*kl+ki;
        extf = sum(sum( (reshape(ppw(1:oi-1),oi-1,1)*ks).^powTime) ) + sum( (ppw(oi)*ks(1:ki)).^powTime);
        
        par.k = ks(ki);
        par.N = par.ppw*par.k;
        par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
        par.colltau = par.t(1:par.N);
        
        %% Computate the full solution
        A1 = zeros(par.N); tic; prevToc = toc;
        for i = 1:par.N
%             if (ki ~= 1) && (toc-prevToc > printtoc)
            if (toc-prevToc > printtoc) % Will print wrong expected end when ki = 1 = oi
                prevToc = toc;
                display(['Obst. ' num2str(obstacle) ' at k=' num2str(ks(ki)) ': ' num2str( (i-1)/par.N,'%7.3f'),...
                    '=A1%, now=' datestr(now) ', est. # sec. left for A1=' num2str(toc*(par.N-i+1)/(i-1))  ', est. end ' datestr(expectedEnd - ...
                    (sum(v.timeA(:,1))/extf*(ppw(oi)*ks(ki)).^powTime - toc*par.N/(i-1) )/24/3600)])
            end
            A1(i,:) = collRowQBF(i,par);
        end
        v.timeA(idx,1) = toc;
        b = par.bc(par.k,par.par(par.colltau)); % b is needed for computing the temporary errors when making A2
        c1 = A1\b; % The solution is needed for computing the correlations for ki==1 and for the temporary errors
        
        %% Compute \tilde{A}
        if ki == 1
            % Re-use the full solution, also for computing R.
            A2 = A1;
            v.timeA(idx,2) = v.timeA(idx,1);
        else
%         elseif any(ks(ki) == ksR)
            % Use allBounds from previous entry in ksR to set up A2
            A2 = zeros(par.N); % Could actually use a sparse matrix structure now because its structure is given by R
%             warning off; % Ignore warnings about removing windows not identically one: because just 1 corr above xi.
            tic;
            prevToc = toc;
            % We used a local threshold (max of each row) and not a global threshold (of all of R)
            for i = 1:par.N
                if (toc-prevToc > printtoc)
                    prevToc = toc;
                    display(['Obst. ' num2str(obstacle) ' at k=' num2str(ks(ki)) ': ' num2str(i/par.N,'%7.3f') ...
                        '=A2%, now=' datestr(now) ', est. # sec. left for A2=' num2str(toc*(par.N-i)/i) ', temp. compr. error=' ...
                        num2str(norm(A2(1:(i-1),:)*c1-b(1:(i-1)))/norm(b(1:(i-1))  )) ...
                        ', est. end ' datestr(expectedEnd - (sum(v.timeA(:,2))/extf*(ppw(oi)*ks(ki)).^powTime - toc*par.N/i)/24/3600)]);
                end
                [~, cli] = min(abs(colLow -par.colltau(i)));
                if 0
                    %                 bounds = [];
                    if any(ks(ki) == ksR) || (ks(ki) > ksR(end))
                        [~, cli] = min(abs(colNow -par.colltau(i)));
                        bounds = allB(:,:,cli);
                        %  A2(i,:) = windRow(-i, par, allB(:,:,cli), 0, 1);
                        % A2(i,:) = windRow(i,par,allB(:,:,i),0,1);
                        % elseif (ks(ki) > ksR(end))
                        %   [~, cli] = min(abs(colNow -par.colltau(i)));
                        %   bounds = allB(:,:,cli);
                    else
                        % Interpolate between allBip and allBold: should be ks(ki)^(-1/(r+1)) but we don't know the
                        % order of the stationary point, so we take the most safe option of linear interpolation as
                        % interpolation for all r > 0 will be between O(1/k) and O(k) because it will give the largest window of both.
                        % ounds = allBold + (ks(ki) - ksRold)/(ksRip - ksRold)*allBip;
                        bounds = allBold(:,:,cli); % To get lower errors but less compression?
                        % bounds = allBold(:,:,cli) + (ks(ki) - ksRold)/(ksRip - ksRold)*(allBip(:,:,cli) ...
                        % -allBold(:,:,cli)); %REMOVE WHERE allBip or allBold are zero
                        % A2(i,:) = windRow(-i, par, bounds, 0, 1);
                    end
                    bounds( :, ~any(bounds,1) ) = []; % Delete zero columns as they do not refer to a window, bounds may be empty after this
                    A2(i,:) = windRow(-i, par, bounds, 0, 1);
                elseif 1
                    % From adapAsyComprSingle:
                    tc = par.colltau(i);
                    curThr = par.xi*max(abs(R(cli,:)));
                    I = find(abs(R(cli,:)) >= curThr);
                    ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
                    
                    % bounds = [(tc -cfTg*(1-pDg)/par.k), sigma(I(ini(1:(length(ini)-1) )+1)); (tc + cfTg*(1-pDg)/par.k), sigma(I(ini(2:length(ini))))];
                    bounds = [(tc -0.02), sigma(I(ini(1:(length(ini)-1) )+1)); (tc +0.02), sigma(I(ini(2:length(ini))))]; % aac
                    % bounds = [tc-Tcor, sigma(I(ini(1:(length(ini)-1) )+1)); tc+Tcor, sigma(I(ini(2:length(ini))))]; % id 1 on tc\pm T, aac
                    
                    decay = repmat(0.02,1,size(bounds,2)); % aac
                    % decay = min(Twi, max(pDg*cfTg/par.k,pbd*(bounds(2,:) -bounds(1,:)) ));
                    % decay = max(pDg*cfTg/par.k,pbd*(bounds(2,:) -bounds(1,:)) );
                    % decay(1) = cfTg*pDg/par.k; % Not necessary because pbd*(bounds(2,:) -bounds(1,:)) will be smaller
                    decay(1) = 0.02*16/ks(ki);%oldAar
                    % decay = repmat(cfTg*pDg/par.k,1,size(bounds,2));
                    
                    % decay = repmat(Tcor*percDecay,1,size(bounds,2)); %aac
                    % decay(1) = Tcor*(ks(1)/ks(ki)); % Increasingly smaller window for Green singularity, aac
                    
                    % bounds = [bounds; [1;1]*decay];
                    bounds = [(bounds + [-decay; decay]); [1;1]*decay];%aac
                    
                    A2(i,:) = windRow(i,par,bounds,0,1);%oldAar
                    % A2(i,:) = windRow(i,par,bounds);%aac
                end
            end
            v.timeA(idx,2) = toc;
%         else
            
%             error('interpol version not implemented');
        end
        %% Compute the correlations where needed
        if any(ks(ki) == ksR)
            % Old way of recompression
            tic;
            Tcor = 0.02;
            [R, sigma] = calcCorr(par, A2\b, Tcor, percDecay, [printtoc,expectedEnd-...
                (sum(v.timeA(:,4))/extf*(ppw(oi)*ks(ki)).^powTime)/24/3600], A2);
%             [R, sigma] = calcCorr(par, c2, cfT/par.k, percDecay, [], A2);
            colLow = par.colltau;
            v.timeA(idx,4) = toc;
        end
        if any(ks(ki) == ksR) && 0
%             if ki > 2
            if ks(ki) > ksR(2)
                allBold = allB;
                colLow = colNow;
                ksRold = ksRip;
            end
            c2 = A2\b;
            tic
            % We will use a local threshold (max of each row) and not a global threshold (of all of R) to compute the bounds
%             allB = calcCorr(par, c2, cfT/ks(ki), percDecay, [printtoc,expectedEnd-...
%                 (sum(v.timeA(:,4))/extf*(ppw(oi)*ks(ki)).^powTime)/24/3600], A2, 'colB', [cfTg, pDg]);
            allB = calcCorr(par, c2, cfT/ks(ki), percDecay, [printtoc,expectedEnd-...
                (sum(v.timeA(:,4))/extf*(ppw(oi)*ks(ki)).^powTime)/24/3600], A2, 'rowB', [cfTg, pDg]);
%             allB = calcCorr(par, c1, cfT/ks(ki), percDecay, [printtoc,expectedEnd-...%Try whether this reduces errors
%                 (sum(v.timeA(:,4))/extf*(ppw(oi)*ks(ki)).^powTime)/24/3600], A1, 'colB', [cfTg, pDg]);
%             colLow = par.colltau; % Save these for higher freqencies.
            
            colNow = par.colltau; % Save these for higher freqencies.
            
            if ki == 1
                allBold = allB;
                colLow = par.colltau; % Save these for higher freqencies.
                ksRold = ks(ki);
            else
                % Find corresponding closest window bounds to allBold
                allBip = 0*allBold; % allB interpolated to allBold
                for i = 1:0%size(allBold,3)
                    [~, cli] = min(abs(colLow(i) -colNow));
                    nzIx = find(allB(1,:,cli));
                    for bi = 1:size(allBold,2)
                        if(allBold(1,bi,i) ~= 0)
%                             maxD = (allBold(2,bi,i) - allBold(1,bi,i))/3; % Bounds should not shift more than this
                            maxD = (allBold(2,bi,i) - allBold(1,bi,i))/2; % Bounds should not shift more than this
%                             [~, ix] = min(abs(allB(1,:,cli) - allBold(1,bi,i)));
%                             [~, ix] = min(abs(allB(1,nzIx,cli) - allBold(1,bi,i)));
%                             mtr = repmat(abs(allB(1,nzIx,cli) - allBold(1,bi,i)), 3,1) + repmat([-1; 0; 1], 1, length(nzIx));
                            mtr = abs(repmat(allB(1,nzIx,cli) - allBold(1,bi,i), 3,1) + repmat([-1; 0; 1], 1, length(nzIx)) );
                            [mix, ixp] = min(mtr, [], 2);
                            [mi, per] = min(mix);
                            ix = ixp(per);
%                             while allB(1,ix,cli) == 0
%                                 [~, ix] = min(abs(allB(1,:,cli) - allBold(1,bi,i))) ?? ;
%                             end
                            if allB(1,ix,cli) +per -2 < allBold(1,bi,i)
                                disp([num2str(ks(ki)) '=k: window is expanding to the left from ' num2str(allBold(1,bi,i))...
                                    ' to ' num2str(allB(1,ix,cli) +per -2)]);
                            end
%                             if (allB(1,ix,cli) +per -2 < allBold(1,bi,i)) || (allB(1,ix,cli) +per -2 > allBold(2,bi,i) )
%                                 continue; % window does not exist any more at this freq
%                             elseif (allB(1,ix,cli) == 0) || any(allBip(1,:,i) == allB(1,ix,cli) +per -2) || (mi > maxD) 
                            if (allB(1,ix,cli) == 0) || any(allBip(1,:,i) == allB(1,ix,cli) +per -2) || (mi > maxD) ...
                                    || (allB(1,ix,cli) +per -2 > allBold(2,bi,i) )
                                error('Wrong compression of bounds')
                            end
                            allBip(1,bi,i) = allB(1,ix,cli) +per -2;
                            allBip(3,bi,i) = allB(3,ix,cli);
                            
%                             [~, ix] = min(abs(allB(2,nzIx,cli) - allBold(2,bi,i)));
                            mtr = abs(repmat(allB(2,nzIx,cli) - allBold(2,bi,i), 3,1) + repmat([-1; 0; 1], 1, length(nzIx)) );
                            [mix, ixp] = min(mtr, [], 2);
%                             [~, per] = min(mix);
                            [mi, per] = min(mix);
                            ix = ixp(per);
%                             allBip(2,bi,i) = allB(2,ix,cli);
                            if allB(2,ix,cli) + per - 2 > allBold(2,bi,i)
                                 disp([num2str(ks(ki)) '=k: window is expanding to the right from ' num2str(allBold(2,bi,i))...
                                    ' to ' num2str(allB(2,ix,cli) + per - 2)]);
                            end
                            if (allB(2,ix,cli) == 0) || any(allBip(2,:,i) == allB(2,ix,cli) +per -2) || (mi > maxD) || ...
                                    (allB(2,ix,cli) + per - 2 < allBold(1,bi,i) )
%                                     (allB(2,ix,cli) + per - 2 < allBold(1,bi,i) ) || (allB(2,ix,cli) + per - 2 > allBold(2,bi,i) )
                                error('Wrong compression of bounds')
                            end
                            allBip(2,bi,i) = allB(2,ix,cli) + per - 2;
                            allBip(4,bi,i) = allB(4,ix,cli);
                        end
                    end
                end
                ksRip = ks(ki);
            end
            v.timeA(idx,4) = toc;
        end
        
        %% Validate
        v = validate(A1,A2,par,v,idx);
%         save('recomprSingle.mat','-regexp','^(?!(A1|A2|R)$).') % Don't save results when plotting sparsity pattern
        v.timeA(idx,3) = (now-startk)*24*3600;
        
        expectedEnd = startAll + (now - startAll)*sum(sum( (ppw*ks).^powTime ) )/extf;
        display(['Obstacle ' num2str(obstacle) ' at k = ' num2str(ks(ki)) ' took ' num2str(v.timeA(idx,3) ) ...
            ' sec and now is ' datestr(now) ', expected end ' datestr(expectedEnd)]);
    end
end

return
%% Validation
fs = 20;
figure;
if exist('R','var')
    % figure; pcolor(log(abs(R)) ); shading interp; set(gca,'YDir','rev'); hh = colorbar(); 
    pcolor(log(max(abs(R),20*par.xi)) ); shading interp; set(gca,'YDir','rev'); hh = colorbar();
    ylabel(hh,'$\log(\max(|R_{m,n}|,20\xi))$', 'interpreter','latex', 'FontSize',fs);
else
    % Plot allBounds
    for bi = 1:size(allB,2)
        plot(squeeze(allB(1,bi,:)), 'b.'); hold on;
        plot(squeeze(allB(2,bi,:)), 'r.');
    end
end
% tb = allBold; figure; for bi = 1:size(tb,2), plot(squeeze(tb(1,bi,:)), 'b.'); hold on; plot(squeeze(tb(2,bi,:)), 'r.');    end
% plotVal(v, 0, {'Circle', 'Ellipse', 'Near-inclusion','Nearly convex', 'Nonconvex polygon'});

return
%% Make plot of the percentages combined with multiple scattering obstacles
vsing = v;
load recomprMultiple.mat
fss = 'Fontsize'; fs = 22;
lws = 'LineWidth'; lw = 5;
mss = 'MarkerSize'; ms = 15;
l = {'+', 'x', 'o', 'v', '*', 'h', 'd'}; 
ll = {'-','--',':', '-.'};
c = 'brgkcmy';

sl = length(vsing.ks);
figure; 
loglog(vsing.ks, vsing.perc(1:sl,2), 'b-', lws, lw); 
hold on;
loglog(vsing.ks, vsing.perc(sl+1:2*sl,2), 'r--', lws, lw); 
loglog(vsing.ks, vsing.perc(2*sl+1:3*sl,2), 'g-.', lws, lw); 
loglog(vsing.ks, vsing.perc(3*sl+1:4*sl,2), 'k:', lws, lw); 
loglog(vsing.ks, vsing.perc(4*sl+1:5*sl,2), 'c--', lws, lw); 

loglog(v.ks, v.perc(1:kl,2), 'm:', lws, lw); 
loglog(v.ks, v.perc(kl+1:2*kl,2), 'b-.', lws, lw); 
loglog(v.ks, v.perc(2*kl+1:3*kl,2), 'k-', lws, lw); 
loglog([vsing.ks(2), vsing.ks(end)], min(min(vsing.perc))*[vsing.ks(2)/vsing.ks(end), 1].^(-1/2), 'r-', lws, lw); 
xlabel('k',fss,fs); 
set(gca,fss,fs);
legend({'Circle', 'Ellipse', 'Near-inclusion','Nearly convex', 'Nonconvex polygon', ...
    '2 Circles','Circle and 2 ellipses', 'Near-inclusion and circle','$\mathcal{O}(k^{-1/2})$'}, 'interpreter', 'latex', fss, 19);


%% Make a plot of the condition numbers
figure; 
loglog(vsing.ks, vsing.conds(2*sl+1:3*sl,2), 'g-.', lws, lw);
hold on;
loglog(vsing.ks, vsing.conds(2*sl+1:3*sl,1), 'g', lws, lw); 
loglog(vsing.ks, vsing.conds(4*sl+1:5*sl,2), 'c--', lws, lw);
loglog(vsing.ks, vsing.conds(4*sl+1:5*sl,1), 'c', lws, lw); 
loglog(v.ks, v.conds(1:kl,2), 'm:', lws, lw);
loglog(v.ks, v.conds(1:kl,1), 'm', lws, lw);
loglog([vsing.ks(1), vsing.ks(end)], max(max(v.conds))*[vsing.ks(1)/vsing.ks(end), 1].^(3/2), 'r-', lws, lw); 
xlabel('k',fss,fs); 
ylabel('Condition number',fss,fs); 
set(gca,fss,fs);
legend({'Near-inclusion $\tilde{A}$','Near-inclusion $A$','Nonconvex polygon $\tilde{A}$', 'Nonconvex polygon $A$', ...
    'Two circles $\tilde{A}$', 'Two circles $A$','$\mathcal{O}(k^{3/2})$'}, 'interpreter', 'latex', fss, 19);

%% Make a plot of the number of GMRES iterations
figure
loglog(vsing.ks, vsing.nbIter(2*sl+1:3*sl,2), 'g-.', lws, lw);
hold on;
loglog(vsing.ks, vsing.nbIter(2*sl+1:3*sl,1), 'g', lws, lw); 
loglog(vsing.ks, vsing.nbIter(4*sl+1:5*sl,2), 'c--', lws, lw);
loglog(vsing.ks, vsing.nbIter(4*sl+1:5*sl,1), 'c', lws, lw); 
loglog(v.ks, v.nbIter(1:kl,2), 'm:', lws, lw);
loglog(v.ks, v.nbIter(1:kl,1), 'm', lws, lw);
loglog([vsing.ks(1), vsing.ks(end)], vsing.nbIter(3*sl,2)*[vsing.ks(1)/vsing.ks(end), 1].^(1/4), 'r-', lws, lw); 
xlabel('k',fss,fs); 
ylabel('# GMRES iterations',fss,fs); 
set(gca,fss,fs);
legend({'Near-inclusion $\tilde{A}$','Near-inclusion $A$','Nonconvex polygon $\tilde{A}$', 'Nonconvex polygon $A$', ...
    'Two circles $\tilde{A}$', 'Two circles $A$','$\mathcal{O}(k^{1/4})$'}, 'interpreter', 'latex', fss, 19);
