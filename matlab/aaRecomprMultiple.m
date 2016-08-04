% Adaptive Asymptotic Recompression for multiple scattering obstacles: recompute correlations at each wavenumber in an automatic and 
% efficient way to get decreasing window function supports for general geometries.
%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

percDecay = 1;
Tcor = 0.02; 
fr = 1.5; % Fraction of N for number of columns of R

avm = 100; % Number of random taus to average BC over
taus = rand(avm,1);

ks = 2.^(4:10);
obsts = 5:7;
kl = length(ks);
maxob = length(obsts);
printtoc = 300;
mti = 2;
v = struct('conds', zeros(maxob*kl,2), 'mti', mti, 'avm', avm, 'taus', taus, 'errBCavm', zeros(maxob*kl,2+mti),...
    'perc', zeros(maxob*kl,2), 'errSol', zeros(maxob*kl,2+mti), 'compresErr', zeros(maxob*kl,2),  'errInt', zeros(maxob*kl,2+mti), ...
    'timeSol', zeros(maxob*kl,2+mti), 'nbIter', zeros(maxob*kl,mti), 'timeA', zeros(maxob*kl,4), 'ks', ks);
ppw = zeros(maxob,1); % Total N divided by k, for the time estimates
for oi = 1:length(obsts)
    par = getObst(obsts(oi));
    for obst = 1:length(par.obsts)
        ppw(oi) = ppw(oi) + par.obsts(obst).ppw;
    end
end
startAll = now;
powTime = 1.5;
for oi = 1:length(obsts)
	obstacle = obsts(oi);
	par = getObst(obstacle);
	for ki = 1:kl
 	       startk = now;
		idx = (oi-1)*kl+ki;
        	TcorN = Tcor*(ks(1)/ks(ki));
		par.k = ks(ki);
		par.N = 0;
		par.r = zeros(2,length(par.obsts)); % ranges
	        for obst = 1:length(par.obsts)
			par.obsts(obst).k = par.k;
			par.obsts(obst).N = par.obsts(obst).ppw*par.k;
			par.r(1,obst) = par.N+1;
			par.N = par.N + par.obsts(obst).N;
			par.r(2,obst) = par.N;
			
			par.obsts(obst).t = linspace(0,1,par.obsts(obst).N+1); % The knots of the periodic spline;
			par.obsts(obst).colltau = par.obsts(obst).t(1:par.obsts(obst).N);
        	end
		
		%% Computating full solution
		A1 = zeros(par.N); tic; prevToc = toc;
		for i = 1:par.N
			if (ki ~= 1) && (toc-prevToc > printtoc)
				prevToc = toc;
				display(['Obst. ' num2str(obstacle) ' at k=' num2str(ks(ki)) ': ' num2str( (i-1)/par.N,'%7.3f'),...
                    '=A1%, now=' datestr(now) ', est. # sec left for A1=' num2str(toc*(par.N-i+1)/(i-1)) ', est. end ' datestr(expectedEnd-...
                    (sum(v.timeA(:,1))/extf*(ppw(oi)*ks(ki)).^powTime - toc*par.N/(i-1) )/24/3600)])
			end
			obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
			collx = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) );
			for obst = 1:length(par.obsts)
				if obst == obsin
					A1(i,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst));
				else
					A1(i,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collx);
				end
			end
		end
		v.timeA(idx,1) = toc;
		% The solution is needed for computing the correlations for ki==1 and for the temporary errors
		b = zeros(par.N,1);
		for obst = 1:length(par.obsts)
			b(par.r(1,obst):par.r(2,obst)) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).colltau));
		end
		c1 = A1\b;
        
        if ki == 1
            % The full solution is already computed: re-use it, also for computing the correlations
            A2 = A1;
            v.timeA(idx,2) = v.timeA(idx,1);
            if oi == 1, expectedEnd = now; extf = 0; end
        else
            %% Computing A2
            curThr = par.xi*max(max(abs(R)));
            A2 = zeros(par.N);
            tic;
            prevToc = toc;
            for i = 1:par.N
                obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
                tc = par.obsts(obsin).colltau(i-par.r(1,obsin)+1);
                collx = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) );
                [~, cli] = min(abs(colLow{obsin}-tc));
                cli = cli + rlow(1,obsin)-1;
                curThr = par.xi*max(abs(R(cli,:)));
                if (toc-prevToc > printtoc)
                    prevToc = toc;
                    display(['Obst ' num2str(obstacle) ' at k=' num2str(ks(ki)) ': ', num2str(i/par.N,'%7.3f') '=A2%, now=' datestr(now)...
                        ', est. # sec. left for A2=' num2str(toc*(par.N-i)/i) ', tmp comprErr=' num2str(norm(A2(1:(i-1),:)*c1-b(1:(i-1)))/...
                        norm(b(1:(i-1))  )) ', est. end ' datestr(expectedEnd-(sum(v.timeA(:,2))/extf*(ppw(oi)*ks(ki)).^powTime - ...
                        toc*par.N/i )/24/3600)]);
                end
                
                for obst = 1:length(par.obsts)
                    rowt = R(cli,:);
                    rowt(1:(obbounds(1,obst)-1)) = 0;
                    rowt((obbounds(2,obst)+1):end) = 0;
                    I = find(abs(rowt) >= curThr);
                    ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
                    if (obst == obsin) && (numel(I) > 1)
                        bounds = [tc-Tcor, sigma(I(ini(1:(length(ini)-1) )+1)); tc+Tcor, sigma(I(ini(2:length(ini))))];
                        decay = repmat(Tcor*percDecay,1,size(bounds,2));
                        decay(1) = TcorN;
                        bounds = [(bounds + [-decay; decay]); [1;1]*decay];
                        A2(i,par.r(1,obsin):par.r(2,obsin)) = windRow(i-par.r(1,obsin)+1,par.obsts(obsin),bounds,0,1);
                    elseif (numel(I) > 1)
                        bounds = [sigma(I(ini(1:(length(ini)-1) )+1)); sigma(I(ini(2:length(ini))))];
                        decay = repmat(Tcor*percDecay,1,size(bounds,2));
                        decay(1) = TcorN;
                        bounds = [(bounds + [-decay; decay]); [1;1]*decay];
                        A2(i,par.r(1,obst):par.r(2,obst)) = windRow('error', par.obsts(obst), bounds,0,1,collx);
                    end
                end
            end
            v.timeA(idx,2) = toc;
        end
        
        %% Computing correlations
        if (ki ~= kl) % Compute the R to be used for the next ki
            tic;
            [R, sigma,obbounds] = calcCorr(par, A2\b, Tcor, percDecay, [printtoc,expectedEnd-...
                (sum(v.timeA(:,4))/extf*(ppw(oi)*ks(ki)).^powTime)/24/3600],A2);
            v.timeA(idx,4) = toc;
            colLow = cell(length(par.obsts),1);
            rlow = par.r;
            for obst = 1:length(par.obsts)
                colLow{obst} = par.obsts(obst).colltau;
            end
        end
		v = validate(A1,A2,par,v,idx);
        save('recomprMultiple.mat','-regexp','^(?!(A1|A2|R)$).')
        v.timeA(idx,3) = (now-startk)*24*3600; % Also has the time for validation
        
        extf = sum(sum( (reshape(ppw(1:oi-1),oi-1,1)*ks).^powTime) ) + sum( (ppw(oi)*ks(1:ki)).^powTime);
        expectedEnd = startAll + (now - startAll)*sum(sum( (ppw*ks).^powTime ) )/extf;
        display(['Obstacle ' num2str(obstacle) ' at k = ' num2str(ks(ki)) ' took ' num2str(v.timeA(idx,3) ) ...
            ' sec and now is ' datestr(now) ', expected end ' datestr(expectedEnd)]);
	end
end

%% Compression percentages areplotted in aaRecomprSingle, print timing table here
plotVal(v,0,{'2 Circles','Circle and 2 ellipses', 'Near-inclusion and circle'});
