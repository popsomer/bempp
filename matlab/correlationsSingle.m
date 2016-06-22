% Compute correlations to automatically determine where to put our window functions for general geometries.
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

ks = 2.^(7:8);%10);

% printtoc = 600; % Print progress info every 5s - > for recompression
obsts = [8 3];
mti = 2;

kl = length(ks);
nbOb = length(obsts);

avm = 100; % Number of random taus to average BC over
v = struct('conds', zeros(nbOb*kl,2), 'mti', mti, 'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(nbOb*kl,2+mti),...
	'perc', zeros(nbOb*kl,2), 'errSol', zeros(nbOb*kl,2+mti), 'compresErr', zeros(nbOb*kl,2),  'errInt', zeros(nbOb*kl,2+mti), ...
	'timeSol', zeros(nbOb*kl,2+mti), 'nbIter', zeros(nbOb*kl,mti), 'timeA', zeros(nbOb*kl,4), 'ks', ks);

%% Computations
for oi = 1:length(obsts)
    obstacle = obsts(oi);
	par = getObst(obstacle);
    start = now;
    
	for ki = 1:kl
		idx = (oi-1)*kl+ki;
		par.k = ks(ki);
		par.N = par.ppw*par.k;
		par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
		par.colltau = par.t(1:par.N);
		
		%% Computating full solution
		A1 = zeros(par.N); tic;
		for i = 1:par.N
		% 	parfor i=1:par.N % Instead of a sequential loop
			A1(i,:) = collRowQBF(i,par);
		end
		v.timeA(idx,1) = toc;
        b = par.bc(par.k,par.par(par.colltau)); % b is also needed for temporary errors when making A2.
        c1 = A1\b; % The solution is needed for computing the correlations.
        
        %% Computing correlations
        if ki == 1 % Re-use the rois/roit that is computed here at higher frequencies.
            tic
            [rois, roit] = calcCorr(par, c1, Tcor, percDecay, [inf,0]); % Don't use A1, but the integral.
            v.timeA(idx,4) = toc;
            colLow = par.colltau; % Save these for higher freqencies.
        end
        
		%% Compute A2
		A2 = zeros(par.N); % We could actually use a sparse matrix structure now because its structure is given by rois.
		curThr = par.xi*max(max(abs(rois))); % With a local threshold, the loops above and below can be joined but less efficient
		tic;
		for i = 1:par.N
			tc = par.colltau(i);
			[~, cli] = min(abs(colLow-tc)); % Index i in rois closest to tc
			if thrType == 'l', curThr = par.xi*max(abs(rois(cli,:))); end
			if corrDist
				toAdd = Cinf(roit,tc-Tcor, tc-Tcor/5, tc+Tcor/5, tc+Tcor,1); % Add to the diagonal to ensure the Green singularity is included.
				A2(i,:) = windRow(i, par, [], [], 0, 1, abs(rois(i,:)/curThr) +toAdd, corrDist, roit);
			else % Physical distance
                I = find(abs(rois(cli,:)) >= curThr);
                ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
                bounds = [tc-Tcor, roit(I(ini(1:(length(ini)-1) )+1)); tc+Tcor, roit(I(ini(2:length(ini))))]; % identically 1 on tc \pm T
                decay = repmat(Tcor*percDecay,1,size(bounds,2));
                bounds = [(bounds + [-decay; decay]); [1;1]*decay];
                A2(i,:) = windRow(i,par,bounds);
			end
		end % loop over row-indices i
		v.timeA(idx,2) = toc;
		
		v = validate(A1,A2,par,v,idx);
% 		save('correlationsSingle.mat','-regexp','^(?!(A1|A2|rois)$).')
		save('correlationsSingle.mat','v')
        
        display([num2str(oi) ' = oi, ki = ' num2str(ki) ', now is ' datestr(now) ', expected end ' datestr(start + ...
            (now-start)*sum(ks.^2)/sum(ks(1:ki).^2)*( length(obsts)-oi + 1) )  ]);
	end
end

%% Show the results
figure; surf(abs(rois), 'EdgeColor','none'); xlabel('n'); ylabel('m'); 

fs = 20;
figure; 
pcolor(min(abs(rois),0.1) );
hold on
shading interp;
xlabel('n'); 
ylabel('m'); 
set(gca,'FontSize',fs);

hh = colorbar(); 
ylabel(hh,'min$(|r_{m,n}|,0.1)$', 'interpreter','latex', 'FontSize',fs); 
set(gca,'YDir','rev')

plotVal(v,0,{'Polygon', 'Near-inclusion'})


