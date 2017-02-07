% Investigate the dependence of errors on the threshold.
%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

percDecay = 1; % Percentage of the window for the C-inf decay: 0 means a block window and 1 means not identically one on any interval.
thrType = 'l'; % Or 'g' for global threshold
Tcor = 0.02;

% If corrDist is nonzero, use the correlations iso the physical distance to determine the window. The value gives
% the percentage below the threshold from where the window is identically zero.
corrDist = 0;
ks = 2^8;
printtoc = 10;
mti = 0;
% mti = 2;
thrp = logspace(-5,0,20);
tl = length(thrp);

avm = 100; % Number of random taus to average BC over
taus = rand(avm,1);
% v = struct('conds', zeros(tl,2), 'mti', mti, 'avm', avm, 'taus', taus, 'errBCavm', zeros(tl,2+mti),...
v = struct('mti', mti, 'avm', avm, 'taus', taus, 'errBCavm', zeros(tl,2+mti),...
    'nnz', zeros(tl,2), 'perc', zeros(tl,2), 'errSol', zeros(tl,2+mti), ...
    'errBCcol', zeros(tl,2+mti), 'compresErr', zeros(tl,2),  'errInt', zeros(tl,2+mti), ...
    'timeSol', zeros(tl,2+mti), 'nbIter', zeros(tl,mti), 'timeA', zeros(tl,4), 'ks', ks);

% par = getObst(1);
par = getObst(4);
par.k = ks;
par.N = par.ppw*par.k;
par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
par.colltau = par.t(1:par.N);

%% Standard computation
A1 = zeros(par.N); tic; prevToc = toc;
for i = 1:par.N
    if (toc-prevToc > printtoc)
        prevToc = toc;
        display([num2str( (i-1)/par.N,'%7.3f') ' = fraction of full matrix, estimated end on ' datestr(now + ...
            toc*((tl+5)*par.N/(i-1) - 1)/24/3600 ) ]); %' asdf ' num2str(toc*((tl+1)*par.N/(i-1) - 1)) ])
    end
    % 	parfor i=1:par.N % Instead of sequential loop
    A1(i,:) = collRowQBF(i,par);
end
v.timeA(1,1) = toc;

%% Computing correlations
b = par.bc(par.k,par.par(par.colltau));
c1 = A1\b; % Solution is needed for computing the correlations
% The time for constructing rois without using A1 is about 3.7 times higher than constructing A1.
tic;
[rois, roit] = calcCorr(par, c1, Tcor, percDecay, [printtoc, now + v.timeA(1,1)*(tl+4)/24/3600] ); 
v.timeA(1,4) = toc;

%% Looping over all threshold percentages
start = now;
for ti = 1:tl
    idx = ti;
    %% Compute A2
    A2 = zeros(par.N);
    curThr = thrp(ti)*max(max(abs(rois)));
    tic;
    prevToc = toc;
    warning off; % Avoid warnings of ill-conditioned systems when the threshold percentage is too low.
    for i = 1:par.N
        if (toc-prevToc > printtoc)
            prevToc = toc;
            display([num2str(thrp(ti)) ' = \xi, ' num2str((i-1)/par.N,'%7.3f') ' = (i-1)/N, ||A c -b||/||b|| up to row i equals ' ...
                num2str(norm(A2(1:(i-1),:)*c1-b(1:(i-1)))/norm(b(1:(i-1))  )) ', now is ' datestr(now) ', estimated end on ' datestr(start + ...
                (now-start)*par.N*tl/(par.N*(ti-1)+i) ) ])
            
        end
        tc = par.colltau(i);
        if thrType == 'l', curThr = thrp(ti)*max(abs(rois(i,:))); end
        if corrDist
            toAdd = Cinf(roit,tc-Tcor, tc-Tcor/5, tc+Tcor/5, tc+Tcor,1); % Add to the diagonal to ensure the Green singularity is included.
            A2(i,:) = windRow(i, par, [], [], abs(rois(i,:)/curThr) +toAdd, corrDist, roit);
        else % Physical distance
            I = find(abs(rois(i,:)) >= curThr);
            ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
            bounds = [tc-Tcor, roit(I(ini(1:(length(ini)-1) )+1)); tc+Tcor, roit(I(ini(2:length(ini))))]; % id 1 on tc\pm T
            decay = repmat(Tcor*percDecay,1,size(bounds,2));
            bounds = [(bounds + [-decay; decay]); [1;1]*decay];
            A2(i,:) = windRow(i,par,bounds);
        end
    end % loop over row-indices i
    v.timeA(idx,2) = toc;
    warning on;
    
    v = validate(A1,A2,par,v,idx);
    save('threshold.mat','-regexp','^(?!(A1|A2|rois)$).')
    display([num2str(thrp(ti)) ' = \xi ended on ' datestr(now) ', expected end ' datestr(start + (now-start)*tl/ti)]);
end
v

%% Make the figure
fss = 'Fontsize'; fs = 22;
lws = 'LineWidth'; lw = 5;
mss = 'MarkerSize'; ms = 15;
l = {'v', '+', 'o', 'x', '*', 'h', 'd'};
ll = {'-','--',':', '-.'};
c = 'bgrkcmy';

figure;
loglog(thrp, v.perc(:,2), 'r-', lws, lw);
hold on
loglog(thrp, v.errSol(:,2)+eps, 'b-', lws,lw)
loglog(thrp, v.compresErr(:,2), 'k-.', lws,lw);
loglog(thrp, v.errBCavm(:,2), 'r--', lws,lw);
loglog(thrp, v.errBCavm(:,1), 'g:', lws,lw);
loglog(thrp, thrp, 'b:', lws, lw);
loglog(thrp, v.errInt(:,2), 'm-', lws,lw);
loglog(thrp, v.errInt(:,1), 'c:', lws,lw);

% legend({'$||c-\tilde{c}||/||c||$', 'err BC $\tilde{c}$', 'err BC $c$', '$||\tilde{A}c-b||/||b||$', '$\xi$'...
%     'int field $\tilde{c}$', 'int field $c$', '\% nnz'}, 'interpreter','latex', fss, fs)
legend({'\% nonzeros', '$||c-\tilde{c}||/||c||$', '$||\tilde{A}c-b||/||b||$', 'Res. i.e. $\tilde{c}$', 'Res. i.e. $c$', '$\xi$'...
    'Int. field $\tilde{c}$', 'Int. field $c$'}, 'interpreter','latex', fss, fs)
xlabel('\xi')
ylabel('Error')
set(gca,fss,fs);
