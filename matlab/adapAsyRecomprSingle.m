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
pDg = 0.8; cfTg = 1.2;

ksR = 2.^(4:7); % The wavenumbers at which to recompute R using A2
% ksR = 2.^(4); % The wavenumbers at which to recompute R using A2
% ksR = 2^9;
dk = 16; % Distance between consecutive wavenumbers
% obsts = [1 2 3 4 8];
obsts = 4; %3;

% ks is a list of all wavenumbers and must start with ksR(1:2), then the wavenumbers where we apply
% interpolation between those, then ksR(3) for interpolation of R at ksR(2) and ks(3) and so on
% ks = [ksR, (ksR(1)+dk):dk:(ksR(2)-dk)];
ks = ksR;
kl = length(ks);
maxob = length(obsts);
mti = 2; % Number of tests using iterative solvers, can be zero
avm = 100; % Number of random taus to average the residue of the integral equation over
v = struct('conds', zeros(maxob*kl,2), 'mti', mti, 'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(maxob*kl,2+mti),...
    'perc', zeros(maxob*kl,2), 'errSol', zeros(maxob*kl,2+mti), 'compresErr', zeros(maxob*kl,2),  'errInt', zeros(maxob*kl,2+mti), ...
    'timeSol', zeros(maxob*kl,2+mti), 'nbIter', zeros(maxob*kl,mti), 'timeA', zeros(maxob*kl,4), 'ks', ks);

printtoc = 3; % Print info after this number of seconds
powTime = 1.5; % The estimated complexity for printing the estimated remaining time
ppw = zeros(maxob,1); % Number of points per wavelength to estimate the remaining time
for oi = 1:length(obsts)
    par = getObst(obsts(oi));
    ppw(oi) = par.ppw;
end


%% Computation
startAll = now;
expectedEnd = startAll;
figure; % For plotting max of corr per row for each ksR
for oi = 1:length(obsts)
    obstacle = obsts(oi);
    par = getObst(obstacle);
    
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
                A2(i,:) = windRow(i,par,allB(:,:,cli),0,1);
%                 A2(i,:) = windRow(i,par,allB(:,:,i),0,1);
            end
        end
        %% Compute the correlations where needed
        if any(ks(ki) == ksR)
            if ki > 2
                allBold = allB;
            end
            c2 = A2\b;
            tic
            % We will use a local threshold (max of each row) and not a global threshold (of all of R) to compute the bounds
%             allB = calcCorr(par, c2, cfT/ks(ki), percDecay, [printtoc,expectedEnd-...
%                 (sum(v.timeA(:,4))/extf*(ppw(oi)*ks(ki)).^powTime)/24/3600], A2, 'colB', [cfTg, pDg]);
            allB = calcCorr(par, c2, cfT/ks(ki), percDecay, [printtoc,expectedEnd-...
                (sum(v.timeA(:,4))/extf*(ppw(oi)*ks(ki)).^powTime)/24/3600], A2, 'rowB', [cfTg, pDg]);
            v.timeA(idx,4) = toc;
            colLow = par.colltau; % Save these for higher freqencies.
            if ki == 1
                allBold = allB;
            end
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
