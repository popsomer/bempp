% Adaptive Asymptotic Recompression for single scattering obstacles: recompute correlations at each wavenumber in an automatic and 
% efficient way to get decreasing window function supports for general geometries.
%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

percDecay = 1; % else R wrong or could define another \xi for A2?
thrType = 'l'; % Or 'g' for global threshold
Tcor = 0.02;
corrDist = 0; % Nonzero then use correlation iso physical distance to determine window, value gives perc below thr where wind=0
ks = 2.^(4:11);
printtoc = 30; 
obsts = [1 2 3 4 8];
mti = 2;
kl = length(ks);
maxob = length(obsts);

avm = 100; % Number of random taus to average BC over
v = struct('conds', zeros(maxob*kl,2), 'mti', mti, 'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(maxob*kl,2+mti),...
    'perc', zeros(maxob*kl,2), 'errSol', zeros(maxob*kl,2+mti), 'compresErr', zeros(maxob*kl,2),  'errInt', zeros(maxob*kl,2+mti), ...
    'timeSol', zeros(maxob*kl,2+mti), 'nbIter', zeros(maxob*kl,mti), 'timeA', zeros(maxob*kl,4), 'ks', ks);
ppw = zeros(maxob,1);
for oi = 1:length(obsts)
    par = getObst(obsts(oi));
    ppw(oi) = par.ppw;
end

startAll = now;
for oi = 1:length(obsts)
    obstacle = obsts(oi);
    par = getObst(obstacle);
    for ki = 1:kl
        startk = now;
        idx = (oi-1)*kl+ki;
        par.k = ks(ki);
        par.N = par.ppw*par.k;
        par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
        par.colltau = par.t(1:par.N);
        TcorN = Tcor*(ks(1)/ks(ki));
        
        %% Computate full solution
        A1 = zeros(par.N); tic; prevToc = toc;
        for i = 1:par.N
            if (ki ~= 1) && (toc-prevToc > printtoc)
                prevToc = toc;
                display(['Obst. ' num2str(obstacle) ' at k=' num2str(ks(ki)) ': ' num2str( (i-1)/par.N,'%7.3f'),...
                    '=A1%, now=' datestr(now) ', est. # sec. left for A1=' num2str(toc*(par.N-i+1)/(i-1))  ', est. end ' datestr(expectedEnd - ...
                    (sum(v.timeA(:,1))/extf*(ppw(oi)*ks(ki)).^powTime - toc*par.N/(i-1) )/24/3600)])
            end
            A1(i,:) = collRowQBF(i,par);
        end
        v.timeA(idx,1) = toc;
        b = par.bc(par.k,par.par(par.colltau)); % b is needed for temporary errors when making A2
        c1 = A1\b; % The solution is needed for computing the correlations for ki==1 and for the temporary errors
        
        %% Compute correlations
        if ki == 1
            tic
            [R, sigma] = calcCorr(par, c1, Tcor, percDecay, [], A1);
            v.timeA(idx,4) = toc;
            colLow = par.colltau; % Save these for the higher freqency.
            A2 = A1;
            v.timeA(idx,2) = v.timeA(idx,1);
        else
            %% Compute A2
            A2 = zeros(par.N); % Could actually use a sparse matrix structure now because its structure is given by R
            curThr = par.xi*max(max(abs(R))); % For a global threshold.
            warning off; % Ignore warnings about removing windows not identically one: because just 1 corr above xi.
            tic;
            prevToc = toc;
            for i = 1:par.N
                if (toc-prevToc > printtoc)
                    prevToc = toc;
                    display(['Obst. ' num2str(obstacle) ' at k=' num2str(ks(ki)) ': ' num2str(i/par.N,'%7.3f') '=A2%, now=' datestr(now) ...
                        ', est. # sec. left for A2=' num2str(toc*(par.N-i)/i) ', temp. compr. error=' num2str(norm(A2(1:(i-1),:)*c1-b(1:(i-1)))/norm(b(1:(i-1))  ))...
                        ', est. end ' datestr(expectedEnd - (sum(v.timeA(:,2))/extf*(ppw(oi)*ks(ki)).^powTime - toc*par.N/i)/24/3600)]);
                end
                tc = par.colltau(i);
                [~, cli] = min(abs(colLow-tc)); % Index i in R closest to tc
                if thrType == 'l', curThr = par.xi*max(abs(R(cli,:))); end
                I = find(abs(R(cli,:)) >= curThr);
                ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
                bounds = [tc-Tcor, sigma(I(ini(1:(length(ini)-1) )+1)); tc+Tcor, sigma(I(ini(2:length(ini))))]; % id 1 on tc\pm T
                decay = repmat(Tcor*percDecay,1,size(bounds,2));
                decay(1) = TcorN; % Increasingly smaller window for Green singularity
                bounds = [(bounds + [-decay; decay]); [1;1]*decay];
                A2(i,:) = windRow(i,par,bounds);
            end % loop over row-indices i
            v.timeA(idx,2) = toc;
            warning on;
        end
        
        %% Recompute correlations
        if (ki ~= kl) && (ki ~= 1) % Compute the R to be used for the next ki
            tic
            [R, sigma] = calcCorr(par, A2\b, Tcor, percDecay, [printtoc,expectedEnd-...
                (sum(v.timeA(:,4))/extf*(ppw(oi)*ks(ki)).^powTime)/24/3600], A2);
            v.timeA(idx,4) = toc;
            colLow = par.colltau; % Save these for higher freqencies.
        end
        v = validate(A1,A2,par,v,idx);
        save('recomprSingle.mat','-regexp','^(?!(A1|A2|R)$).')
        v.timeA(idx,3) = (now-startk)*24*3600;
        
        powTime = 1.5;
        extf = sum(sum( (reshape(ppw(1:oi-1),oi-1,1)*ks).^powTime) ) + sum( (ppw(oi)*ks(1:ki)).^powTime);
        expectedEnd = startAll + (now - startAll)*sum(sum( (ppw*ks).^powTime ) )/extf;
        display(['Obstacle ' num2str(obstacle) ' at k = ' num2str(ks(ki)) ' took ' num2str(v.timeA(idx,3) ) ...
            ' sec and now is ' datestr(now) ', expected end ' datestr(expectedEnd)]);
    end
end

%% Print a table of the error on the boundary conditions
plotVal(v, 0, {'Circle', 'Ellipse', 'Near-inclusion','Nearly convex', 'Nonconvex polygon'});

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
