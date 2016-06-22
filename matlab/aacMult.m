% Test adaptive asymptotic compression on multiple scattering obstacles
%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

percDecay = 1; % Percentage of the window for the C-inf decay: 0 means a block window and 1 means not identically one on any interval.
thrType = 'l'; % Or 'g' for global threshold
Tcor = 0.02;

avm = 100; % Number of random taus to average BC over

ks =  2^7; %2.^(7:9);
obsts = 5; % Also obstacle 6 for correlations plot
kl = length(ks);
maxob = length(obsts);

mti = 0;
v = struct('conds', zeros(maxob*kl,2), 'mti', mti, 'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(maxob*kl,2+mti),...
    'nnz', zeros(maxob*kl,2), 'perc', zeros(maxob*kl,2), 'errSol', zeros(maxob*kl,2+mti), ...
    'errBCcol', zeros(maxob*kl,2+mti), 'compresErr', zeros(maxob*kl,2),  'errInt', zeros(maxob*kl,2+mti), ...
    'timeSol', zeros(maxob*kl,2+mti), 'nbIter', zeros(maxob*kl,mti), 'timeA', zeros(maxob*kl,4), 'ks', ks);

%% Computations
for oi = 1:length(obsts)
	obstacle = obsts(oi);
	par = getObst(obstacle);
    start = now;
	for ki = 1:kl
        startk = now;
		idx = (oi-1)*kl+ki;
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
		A1 = zeros(par.N); tic;
		for i = 1:par.N
			obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
			collx = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) );
			for obst = 1:length(par.obsts)
				if obst == obsin
					A1(i,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N);
				else
					A1(i,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N,  [], collx);
				end
			end
		end
		v.timeA(idx,1) = toc;
        
		% Solution is needed for computing the correlations
		b = zeros(par.N,1);
		for obst = 1:length(par.obsts)
			b(par.r(1,obst):par.r(2,obst)) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).colltau));
		end
		c1 = A1\b;
        
        %% Computing correlations
        if ki == 1
            tic;
%             [rois, roit,obbounds] = calcCorr(par, c1, Tcor, percDecay, [inf,0]); % Don't use A1, but the integral.
            [rois, roit,obbounds] = calcCorr(par, c1, Tcor, percDecay, [2,now]); % Don't use A1, but the integral.
            v.timeA(idx,4) = toc;
            colLow = cell(length(par.obsts),1);
            rlow = par.r;
            for obst = 1:length(par.obsts)
                colLow{obst} = par.obsts(obst).colltau;
            end
        end
        
		%% Computing A2
        curThr = par.xi*max(abs(rois));
		A2 = zeros(par.N);
		tic;
		prevToc = toc;
		for i = 1:par.N
			obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
            tc = par.obsts(obsin).colltau(i-par.r(1,obsin)+1);
			collx = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) );
            [~, cli] = min(abs(colLow{obsin}-tc));
            cli = cli + rlow(1,obsin)-1;
            if thrType == 'l', curThr = par.xi*max(abs(rois(cli,:))); end
			
			for obst = 1:length(par.obsts)
				rowt = rois(cli,:);
				rowt(1:(obbounds(1,obst)-1)) = 0;
				rowt((obbounds(2,obst)+1):end) = 0;
				I = find(abs(rowt) >= curThr);
				ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
				if (obst == obsin) && (numel(I) > 1)
					bounds = [tc-Tcor, roit(I(ini(1:(length(ini)-1) )+1)); tc+Tcor, roit(I(ini(2:length(ini))))];
                    decay = repmat(Tcor*percDecay,1,size(bounds,2)); 
                    bounds = [(bounds + [-decay; decay]); [1;1]*decay];
					A2(i,par.r(1,obsin):par.r(2,obsin)) = windRow(i-par.r(1,obsin)+1,par.obsts(obsin),bounds);
				elseif (numel(I) > 1)
					bounds = [roit(I(ini(1:(length(ini)-1) )+1)); roit(I(ini(2:length(ini))))];
                    decay = repmat(Tcor*percDecay,1,size(bounds,2)); 
                    bounds = [(bounds + [-decay; decay]); [1;1]*decay];
					A2(i,par.r(1,obst):par.r(2,obst)) = windRow('error', par.obsts(obst), bounds,[],[],collx);
				end
			end
		end
		v.timeA(idx,2) = toc;
		v = validate(A1,A2,par,v,idx);
        v.timeA(idx,3) = (now-startk)*24*3600;
        save('aacMult.mat','v');
        
        display([num2str(oi) ' = oi, ki = ' num2str(ki) ', now is ' datestr(now) ', expected end ' datestr(start + ...
            (now-start)*sum(ks.^2)/sum(ks(1:ki).^2)*( length(obsts)-oi + 1) )  ]);
	end
end


%% Plot spy for article
figure;
fs = 20;
spy(A2);
hold on;
par.obsts(3).colltau(par.N*5/6-par.r(1,3)+2)
h = text(par.N*5/6+2-40, -100, '$\vec{p}$', 'interpreter', 'latex');
set(h,'FontSize',fs, 'color', 'r');
plot(par.N*5/6+2*ones(2,1), [-50; par.N+120], 'r', 'LineWidth', 3);
ylabel('i'); xlabel('j');
set(gca,'FontSize',fs);

%% Plot rois for article
fs = 20;
figure; 
pcolor(min(abs(rois),0.1) );
hold on
shading interp; 
xlabel('n'); 
ylabel('m'); 
set(gca,'FontSize',fs);

h = text(size(rois,2)*5/6+2-40, -100, '$\vec{p}$', 'interpreter', 'latex');
set(h,'FontSize',fs, 'color', 'r');
plot(size(rois,2)*5/6+2*ones(2,1), [-50; size(rois,2)+120], 'r', 'LineWidth', 3);

hh = colorbar(); 
ylabel(hh,'min$(|r_{m,n}|,0.1)$', 'interpreter','latex', 'FontSize',fs); 
set(gca,'YDir','rev')

plotVal(v,0,{'Two circles', 'Three ellipses'})

