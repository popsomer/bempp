% Adaptive Asymptotic Compression for single scattering obstacles: compute correlations to 
% automatically determine where to put our window functions for general geometries.

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

ks = 2^7;
obsts = 3;
mti = 0;

kl = length(ks);
nbOb = length(obsts);

avm = 100; % Number of random taus to average BC over
v = struct('mti', mti, 'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(nbOb*kl,2+mti),...
	'perc', zeros(nbOb*kl,2), 'errSol', zeros(nbOb*kl,2+mti), 'compresErr', zeros(nbOb*kl,2),  'errInt', zeros(nbOb*kl,2+mti), ...
	'timeSol', zeros(nbOb*kl,2+mti), 'nbIter', zeros(nbOb*kl,mti), 'timeA', zeros(nbOb*kl,4), 'ks', ks);

%% Computations
for oi = 1:length(obsts)
	obstacle = obsts(oi);
	par = getObst(obstacle);
	start = now;
    
	for ki = 1:kl-(kl-1)*(obstacle == 3)
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
        b = par.bc(par.k,par.par(par.colltau));
        c1 = A1\b; % The solution is needed for computing the correlations.
        
        %% Computing correlations
        if ki == 1 % Re-use the R that is computed here at higher frequencies.
            tic
            [R, sigma] = calcCorr(par, c1, Tcor, percDecay); % Don't use A1, but the integral.
            v.timeA(idx,4) = toc;
            colLow = par.colltau; % Save these for higher freqencies.
        end
        
        %% Compute A2
        A2 = zeros(par.N); % We could actually use a sparse matrix structure now because its structure is given by R.
		curThr = par.xi*max(max(abs(R)));
		tic;
		for i = 1:par.N
			tc = par.colltau(i);
			[~, cli] = min(abs(colLow-tc)); % Index i in R closest to tc
			if thrType == 'l', curThr = par.xi*max(abs(R(cli,:))); end
			if corrDist
				toAdd = Cinf(sigma,tc-Tcor, tc-Tcor/5, tc+Tcor/5, tc+Tcor,1); % Add to the diagonal to ensure the Green singularity is included.
				A2(i,:) = windRow(i, par, [], [], 0, 1, abs(R(i,:)/curThr) +toAdd, corrDist, sigma);
            else % Physical distance
                I = find(abs(R(cli,:)) >= curThr);
                ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
                bounds = [tc-Tcor, sigma(I(ini(1:(length(ini)-1) )+1)); tc+Tcor, sigma(I(ini(2:length(ini))))]; % identically 1 on tc \pm T
                decay = repmat(Tcor*percDecay,1,size(bounds,2));
                bounds = [(bounds + [-decay; decay]); [1;1]*decay];
                A2(i,:) = windRow(i,par,bounds);
			end
		end % loop over row-indices i
		v.timeA(idx,2) = toc;
		
		v = validate(A1,A2,par,v,idx);
		save('aacSingle.mat','v')
        display([num2str(oi) ' = oi, ki = ' num2str(ki) ', now is ' datestr(now) ', expected end ' datestr(start + ...
            (now-start)*sum(ks.^2)/sum(ks(1:ki).^2)*( length(obsts)-oi + 1) )  ]);
	end
end

%% Show the correlations of the near-inclusion
[R, sigma] = calcCorr(par, c1, Tcor, percDecay, [5, now()], A1); % Approximation of R used to make figure

fs = 20;
figure; 
pcolor(log(max(abs(R),1e-2)) );
hold on
shading interp;
xlabel('n'); 
ylabel('m'); 
set(gca,'FontSize',fs);
lw = 2.7; co = '--r';
h = text(round(size(R,2)*0.35)+2-40, -100, '$\mathbf{q}$', 'interpreter', 'latex');
set(h,'FontSize',fs, 'color', 'r');
plot(round(size(R,2)*0.35)+2*ones(2,1), [-50; size(R,2)+120], co, 'LineWidth', lw);

h = text(round(size(R,2)*0.5)+2-40, -100, '$\mathbf{r}$', 'interpreter', 'latex');
set(h,'FontSize',fs, 'color', 'r');
plot(round(size(R,2)*0.5)+2*ones(2,1), [-50; size(R,2)+120], co, 'LineWidth', lw);

h = text(round(size(R,2)*0.65)+2-40, -100, '$\mathbf{s}$', 'interpreter', 'latex');
set(h,'FontSize',fs, 'color', 'r');
plot(round(size(R,2)*0.65)+2*ones(2,1), [-50; size(R,2)+120], co, 'LineWidth', lw);

hh = colorbar(); 
ylabel(hh,'$\log(\max(|R_{m,n}|,10^{-2}))$', 'interpreter','latex', 'FontSize',fs); 
set(gca,'YDir','rev')

