% Compute correlations to automatically determine where to put our window functions for general geometries.
%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

percDecay = 1; % This is used for computing the correlations and the windows

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
v = struct('conds', zeros(nbOb*kl,2), 'nbGm', mti, 'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(nbOb*kl,2+mti),...
	'perc', zeros(nbOb*kl,2), 'errSol', zeros(nbOb*kl,2+mti), 'compresErr', zeros(nbOb*kl,2),  'errInt', zeros(nbOb*kl,2+mti), ...
	'timeSol', zeros(nbOb*kl,2+mti), 'nbIter', zeros(nbOb*kl,mti), 'timeA', zeros(nbOb*kl,4), 'ks', ks);

%% Computations
expectedEnd = 0;
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
		
		%% Computations
		A1 = zeros(par.N); tic;
		for i = 1:par.N
		% 	parfor i=1:par.N % Instead of sequential loop
			A1(i,:) = collRowQBF(i,1,par.N,par);
		end
		v.timeA(idx,1) = toc;
		
		%% Computing correlations
		b = par.bc(par.k,par.par(par.colltau)); % b is also needed for temporary errors when making A2.
		c1 = A1\b; % The solution is needed for computing the correlations.
        
		if ki == 1 % Re-use the rois/roit that is computed here at higher frequencies.
			rois = zeros(par.N,round(par.N*1.5)); % Compute correlations on more than N points for possible accuracy.
			roit = linspace(0,1,size(rois,2) );
			colLow = par.colltau; % Save these for higher freqencies.
			
			sr = size(rois,2);
            tic;
            % We could compute rois column by column because then the window around roit(roi) would be constant, so that we do not have to
            % recompute the same window indices j1 and j2 for each row. But this would need to recompute the integrand with the 
            % (expensive) Bessel function for each column.
            
            u = [par.t(end-dbf:end-1)-1, par.t, par.t(2:dbf+1)+1]; % the extended knots
            c1ip = deboor(dbf, u, [c1; c1(1:dbf)], roit); % [TODO] Calculate the interpolated solutio
            n
            % Find the intervals of u where the respective roit lies.
            js = [arrayfun(@(roi) find(roi >= u, 1, 'last') - par.dbf - 1, roit), length(u) - 2*par.dbf - 2];
            
            
            for i = 1:par.N
                tc = par.colltau(i);
                integrand = 1i/4.*besselh(0, 1, par.k*sqrt(sum((par.par(tc*ones(size(roit))) ...
                    -par.par(roit) ).^2, 1)) ).*c1ip.*par.gradnorm(roit);
                integrand(isnan(integrand)) = 0; % Temporary fix from taac
                for roi = 1:size(rois,2)
                    [j1, j2, noSplit] = bounds2Ind(roit,roit(roi)-Tcor,roit(roi)+Tcor);
                    if noSplit
                        wi = chi(roit(j1:j2),roit(roi)-Tcor, roit(roi)-(1-percDecay)*Tcor, roit(roi)+Tcor*(1-percDecay),roit(roi)+Tcor,0);
                        rois(i,roi) = sum(wi.*integrand(j1:j2) );
                    else
                        wi1 = chi(roit(j1:sr),roit(roi)-Tcor, roit(roi)-(1-percDecay)*Tcor, roit(roi)+Tcor*(1-percDecay),roit(roi)+Tcor,1);
                        wi2 = chi(roit(1:j2), roit(roi)-Tcor, roit(roi)-(1-percDecay)*Tcor, roit(roi)+Tcor*(1-percDecay),roit(roi)+Tcor,1);
                        rois(i,roi) = sum(wi1.*integrand(j1:sr) ) + sum(wi2.*integrand(1:j2) );
                    end
                end
            end
			v.timeA(idx,4) = toc;
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
				toAdd = Cinf(roit,tc-Tcor, tc-Tcor/5, tc+Tcor/5, tc+Tcor,1);
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
		save('correlationsSingle.mat','-regexp','^(?!(A1|A2|rois)$).')
        
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



