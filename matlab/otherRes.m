% Get additional robustness results.

%% Preliminary test on an integral equation of the second kind
clearvars, close all, format longe
set(0,'DefaultFigureWindowStyle','docked');

par = getObst(1);
% par = getObst(4);
par.secondKind = [1/2, 1]; % This means skpar_a = 1/2, skpar_b = 1: A1 = MM/2 + A1;

par.bc = @(k,x) -1*besselh(0,1, k*sqrt( (x(1,:)' - 0.05).^2 + (x(2,:)' -0.0).^2) ); % Pt source inside circle

Tcor = 0.02;
percDecay = 1;
avm = 100;
printtoc = 3;

par.k = 2^7;
par.N = ceil(par.ppw*par.k);
par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
par.colltau = par.t(1:par.N);

% Computations for the full matrix
A1 = zeros(par.N);
tic;
prevToc = toc;
for i = 1:par.N
    if toc-prevToc > printtoc
        prevToc = toc;
        display([num2str( (i-1)/par.N,'%7.3f') '=A1%, now=' datestr(now) ', est. # sec left for A1=' num2str(toc*(par.N-i+1)/(i-1)) ])
    end
    A1(i,:) = collRowQBF(i,par);
end

MM = zeros(par.N); % Mass matrix
u = [par.t(end-par.dbf:end-1)-1, par.t, par.t(2:par.dbf+1)+1]; % the extended knots
for i=1:par.N
    tc = par.colltau(i);
    for k=1:par.N 
        if par.dbf == 1
            MM(i,k) = MM(i,k) + (tc -u(k))/(u(k+1) -u(k))*((tc >= u(k)) && (tc < u(k+1))) ...
                + (u(k+2) -tc)/(u(k+2) -u(k+1))*((tc >= u(k+1)) && (tc < u(k+2)));
        else
            MM(i,k) = MM(i,k) + bspline_eval(u(k:k+par.dbf+1), par.dbf, tc); % see collRowQBF
        end
		% add the periodic parts for the basis functions at the end
		if k <= par.dbf
            if par.dbf == 1
                MM(i,k) = MM(i,k) + (tc -u(k+par.N))/(u(k+par.N+1) -u(k+par.N))*((tc >= u(k+par.N)) && (tc < u(k+par.N+1))) ...
                    + (u(k+par.N+2) -tc)/(u(k+par.N+2) -u(k+par.N+1))*((tc >= u(k+par.N+1)) && (tc < u(k+par.N+2)));
            else
                MM(i,k) = MM(i,k) + bspline_eval(bfinfo.u(k+par.N:k+par.N+par.dbf+1), par.dbf, tc);
            end
		end
    end
	
end

% Compute the correlations
b = par.bc(par.k,par.par(par.colltau));
c1 = A1\b;
[R, sigma] = calcCorr(par, c1, Tcor, percDecay, [printtoc, now], A1);
colLow = par.colltau;

% Computations for asymptotically compressed matrix
A2 = zeros(par.N);
tic;
prevToc = toc;
for i = 1:par.N
    if toc-prevToc > printtoc
        prevToc = toc;
        display([num2str( (i-1)/par.N,'%7.3f') '=A2%, now=' datestr(now) ', est. # sec left for A2=' num2str(toc*(par.N-i+1)/(i-1))...
            ', tmp comprErr=' num2str(norm(A2(1:(i-1),:)*c1-b(1:(i-1)))/norm(b(1:(i-1))  ))])
    end
    tc = par.colltau(i);
    if 1 % Use correlations
        [~, cli] = min(abs(colLow-tc));
        curThr = par.xi*max(abs(R(cli,:)));
        I = find(abs(R(cli,:)) >= curThr);
        ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
        bounds = [tc-Tcor, sigma(I(ini(1:(length(ini)-1) )+1)); tc+Tcor, sigma(I(ini(2:length(ini))))]; % identically 1 on tc \pm T
        decay = repmat(Tcor*percDecay,1,size(bounds,2));
        bounds = [(bounds + [-decay; decay]); [1;1]*decay];
    else % Capture the SP, specifically for the circle: actually a bit more accurate for same compression
        Tsp = 0.2*(par.k/64)^(-1/2); Tsi = 0.1*(par.k/64)^(-1/2);% So that about 2.7% solerr at k=64
        if (tc < 0.5)
            bounds = [tc-Tsi, 0.5+tc-Tsp; tc+Tsi, 0.5+tc+Tsp; Tsi, Tsp; Tsi, Tsp]; % Singularity and stationary point
        else
            bounds = [tc-0.5-Tsp, tc-Tsi; tc-0.5+Tsp, tc+Tsi; Tsp, Tsi; Tsp, Tsi]; % Singularity and stationary point
        end
        bounds(3:4,:) = bounds(3:4,:)*0.999; %This percentage is for the decay of the windows
    end
    A2(i,:) = windRow(i,par,bounds);
end % loop over row-indices i
c2 = A2\b;

taus = linspace(0+1/avm,1-1/avm,avm);
vSecKind = validate(A1,A2,par,struct('perc', zeros(1,2), 'errSol', zeros(1,2), 'conds', zeros(1,2), 'ks', par.k, 'errInt', zeros(1,2), ...
    'avm', avm, 'taus', taus, 'errBCavm', zeros(1,2) ,'errBCcol', zeros(1,2), 'compresErr', zeros(1,2), 'field', zeros(30) ), 1);
v = struct();
v.field = vSecKind.field;
v.parametr = vSecKind.parametr;
v.xs = vSecKind.xs;
v.ys = vSecKind.ys;
plotVal(v,A2)

vSecKind.errExt = zeros(1,2);
extXs = [(0.8*cos(2*pi*taus)); (0.8*sin(2*pi*taus))]; % Encircles the circle
sgn = +1;
asdf = zeros(3, vSecKind.avm);
for av = 1:vSecKind.avm
    asdf(1,av) = par.bc(par.k, extXs(:,av));
    asdf(2,av) = evalFieldQBF(par, extXs(:,av) ,c1, 0);
    asdf(3,av) = evalFieldQBF(par, extXs(:,av) ,c2, 1);
    vSecKind.errExt(1,1) = vSecKind.errExt(1,1) + abs(par.bc(par.k,extXs(:,av)) - sgn*evalFieldQBF(par,extXs(:,av),c1,0) )/vSecKind.avm;
    vSecKind.errExt(1,2) = vSecKind.errExt(1,2) + abs(par.bc(par.k,extXs(:,av)) - sgn*evalFieldQBF(par,extXs(:,av),c2,1) )/vSecKind.avm;
end
vSecKind.errExt = vSecKind.errExt/mean(abs(asdf(1,:)));
figure; plot(taus, [abs(asdf(2,:)); abs(asdf(2,:) -asdf(1,:)); abs(asdf(3,:) -asdf(1,:))] );
vSecKind.xi = par.xi  % The error on the boundary conditions has not been adjusted to IE of the second kind/Neumann BC


%% Test near resonance frequency of a circle
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

par = getObst(1);
Tcor = 0.02;
percDecay = 1;
avm = 100;
printtoc = 5;

par.k = 2*5.830649350672455e+01; % The first zero of J_{50}(x) divided by the radius of the circle (1/2)
par.N = ceil(par.ppw*par.k);
par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
par.colltau = par.t(1:par.N);

% Full matrix
A1 = zeros(par.N);
tic;
for i = 1:par.N
    A1(i,:) = collRowQBF(i,par);
end

% Correlations
b = par.bc(par.k,par.par(par.colltau));
c1 = A1\b;
[R, sigma] = calcCorr(par, c1, Tcor, percDecay, [], A1);
colLow = par.colltau;

% Asymptotically compressed matrix
A2 = zeros(par.N);
tic;
prevToc = toc;
for i = 1:par.N
    if (toc-prevToc > printtoc)
        prevToc = toc;
        display([num2str(ki,'%8.4f') '=ki, ' num2str(i/par.N,'%7.3f') 'A2%, est. sec left=' num2str(toc*(par.N-i)/i)])
    end
    tc = par.colltau(i);
    [~, cli] = min(abs(colLow-tc));
    curThr = par.xi*max(abs(R(cli,:)));
    I = find(abs(R(cli,:)) >= curThr);
    ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
    bounds = [tc-Tcor, sigma(I(ini(1:(length(ini)-1) )+1)); tc+Tcor, sigma(I(ini(2:length(ini))))]; % identically 1 on tc \pm T
    decay = repmat(Tcor*percDecay,1,size(bounds,2));
    bounds = [(bounds + [-decay; decay]); [1;1]*decay];
    A2(i,:) = windRow(i,par,bounds);
end % loop over row-indices i

% Validation
vAac = validate(A1,A2,par,struct('perc', zeros(1,2), 'errSol', zeros(1,2), 'conds', zeros(1,2), ...
    'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(1,2) ),1)


%% Block window
clearvars, close all, format longe
set(0,'DefaultFigureWindowStyle','docked');

par = getObst(3);
Tcor = 0.02;
percDecay = 1;
avm = 100;

par.k = 128;
par.N = ceil(par.ppw*par.k);
par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
par.colltau = par.t(1:par.N);

% Computations for the full matrix
A1 = zeros(par.N);
tic;
prevToc = toc;
for i = 1:par.N
    A1(i,:) = collRowQBF(i,par);
end

% Compute the correlations
b = par.bc(par.k,par.par(par.colltau));
c1 = A1\b;
[R, sigma] = calcCorr(par, c1, Tcor, percDecay, [], A1);
colLow = par.colltau;

% Computations for standard asymptotically compressed matrix
A2 = zeros(par.N);
tic;
for i = 1:par.N
    tc = par.colltau(i);
    [~, cli] = min(abs(colLow-tc));
    curThr = par.xi*max(abs(R(cli,:)));
    I = find(abs(R(cli,:)) >= curThr);
    ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
    bounds = [tc-Tcor, sigma(I(ini(1:(length(ini)-1) )+1)); tc+Tcor, sigma(I(ini(2:length(ini))))]; % identically 1 on tc \pm T
    decay = repmat(Tcor*percDecay,1,size(bounds,2));
    bounds = [(bounds + [-decay; decay]); [1;1]*decay];
    A2(i,:) = windRow(i,par,bounds);
end % loop over row-indices i

taus = rand(avm,1);
vAac = validate(A1,A2,par,struct('perc', zeros(1,2), 'errSol', zeros(1,2), 'conds', zeros(1,2), ...
    'avm', avm, 'taus', taus, 'errBCavm', zeros(1,2) ),1)

% Block window
A3 = A1;
A3(A2 == 0) = 0;
vBlock = validate(A1,A3,par,struct('perc', zeros(1,2), 'errSol', zeros(1,2), 'conds', zeros(1,2), ...
    'avm', avm, 'taus', taus, 'errBCavm', zeros(1,2) ),1)


%% Continue with cubic basis functions
par.dbf = 3;
A1 = zeros(par.N);
tic;
prevToc = toc;
for i = 1:par.N
    A1(i,:) = collRowQBF(i,par);
end

% Compute the correlations for cubic
b = par.bc(par.k,par.par(par.colltau));
c1 = A1\b;
[R, sigma] = calcCorr(par, c1, Tcor, percDecay, [], A1);
colLow = par.colltau;

% Computations for standard asymptotically compressed matrix for cubic
A2 = zeros(par.N);
tic;
for i = 1:par.N
    tc = par.colltau(i);
    [~, cli] = min(abs(colLow-tc));
    curThr = par.xi*max(abs(R(cli,:)));
    I = find(abs(R(cli,:)) >= curThr);
    ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
    bounds = [tc-Tcor, sigma(I(ini(1:(length(ini)-1) )+1)); tc+Tcor, sigma(I(ini(2:length(ini))))]; % identically 1 on tc \pm T
    decay = repmat(Tcor*percDecay,1,size(bounds,2));
    bounds = [(bounds + [-decay; decay]); [1;1]*decay];
    A2(i,:) = windRow(i,par,bounds);
end 
vCub = validate(A1,A2,par,struct('perc', zeros(1,2), 'errSol', zeros(1,2), 'conds', zeros(1,2), ...
    'avm', avm, 'taus', taus, 'errBCavm', zeros(1,2) ),1)

