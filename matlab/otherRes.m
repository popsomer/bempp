% Use the visibility criterion to determine which regions to include in the windows for the near-inclusion obstacle.

%% Integral equation of the second kind
% clearvars, close all, format longe
set(0,'DefaultFigureWindowStyle','docked');

par = getObst(4);
% par.secondKind = 'second';
par.secondKind = [1/2, 1]; % skpar_a = 1/2, skpar_b = 1: A1 = MM/2 + A1;
par.xi = 0.005;
par.ppw = 14;
par.bc = @(k,x) -1*besselh(0,1, k*sqrt( (x(1,:)' - 0.51).^2 + (x(2,:)' -0.53).^2) ); % Pt source
% par.bc = @(k,x) 1*besselh(0,1, k*sqrt( (x(1,:)' - 0.51).^2 + (x(2,:)' -0.53).^2) ); % Pt source
% dbcx = @(k,x) besselh(1,1, k*sqrt( (x(1,:)' - 0.51).^2 + (x(2,:)' -0.53).^2) ).*...
%     k./2./sqrt( (x(1,:)' - 0.51).^2 + (x(2,:)' -0.53).^2).*2.*(x(1,:)' - 0.51); 
% dbcy = @(k,x) besselh(1,1, k*sqrt( (x(1,:)' - 0.51).^2 + (x(2,:)' -0.53).^2) ).*...
%     k./2./sqrt( (x(1,:)' - 0.51).^2 + (x(2,:)' -0.53).^2).*2.*(x(2,:)' - 0.53); 

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
% skpar_a = 1/2, skpar_b = 1
% A1 = MM/2 + A1;

% Compute the correlations
b = par.bc(par.k,par.par(par.colltau));
% b = -par.bc(par.k,par.par(par.colltau));
% b = -[dbcx(par.k,par.par(par.colltau)), dbcy(par.k,par.par(par.colltau))]*par.normal(par.colltau);
% nor = par.normal(par.colltau);
% b = -(dbcx(par.k,par.par(par.colltau)).*transpose(nor(1,:)) +...
%     dbcy(par.k,par.par(par.colltau)).*transpose(nor(2,:)) );
c1 = A1\b;
[R, sigma] = calcCorr(par, c1, Tcor, percDecay, [printtoc, now], A1);
colLow = par.colltau;

% Computations for standard asymptotically compressed matrix
A2 = A1;%zeros(par.N);
tic;
prevToc = toc;
for i = 1:par.N
    if toc-prevToc > printtoc
        prevToc = toc;
        display([num2str( (i-1)/par.N,'%7.3f') '=A2%, now=' datestr(now) ', est. # sec left for A2=' num2str(toc*(par.N-i+1)/(i-1))...
            ', tmp comprErr=' num2str(norm(A2(1:(i-1),:)*c1-b(1:(i-1)))/norm(b(1:(i-1))  ))])
    end
    tc = par.colltau(i);
    [~, cli] = min(abs(colLow-tc));
    curThr = par.xi*max(abs(R(cli,:)));
    I = find(abs(R(cli,:)) >= curThr);
    ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
    bounds = [tc-Tcor, sigma(I(ini(1:(length(ini)-1) )+1)); tc+Tcor, sigma(I(ini(2:length(ini))))]; % identically 1 on tc \pm T
    decay = repmat(Tcor*percDecay,1,size(bounds,2));
    bounds = [(bounds + [-decay; decay]); [1;1]*decay];
%     A2(i,:) = windRow(i,par,bounds);
    A2(i,:) = windRow(i, par, bounds, 0, 1);
end % loop over row-indices i
% A2 = MM/2 +A2;
c2 = A2\b;

% taus = rand(avm,1);
taus = linspace(0+1/avm,1-1/avm,avm);
% vSecKind = validate(A1,A2,par,struct('perc', zeros(1,2), 'errSol', zeros(1,2), 'conds', zeros(1,2), ...
%     'avm', avm, 'taus', taus, 'errBCavm', zeros(1,2) ),1)
vSecKind = validate(A1,A2,par,struct('perc', zeros(1,2), 'errSol', zeros(1,2), 'conds', zeros(1,2), 'ks', par.k, ...
    'avm', avm, 'taus', taus, 'errBCavm', zeros(1,2) ,'errBCcol', zeros(1,2), 'compresErr', zeros(1,2), 'field', zeros(30) ), 1);
v = struct();
v.field = vSecKind.field;
v.parametr = vSecKind.parametr;
v.xs = vSecKind.xs;
v.ys = vSecKind.ys;
plotVal(v,A2)

vSecKind.errExt = zeros(1,2);
extXs = [(0.5 + 0.6*cos(2*pi*taus)); (0.5 + 0.6*sin(2*pi*taus))]; %Encircles nearly convex obstacle
asdf = zeros(3, vSecKind.avm);
for av = 1:vSecKind.avm
    asdf(1,av) = par.bc(par.k, extXs(:,av));
    asdf(2,av) = evalFieldQBF(par, extXs(:,av) ,c1, 0);
    asdf(3,av) = evalFieldQBF(par, extXs(:,av) ,c2, 1);
    vSecKind.errExt(1,1) = vSecKind.errExt(1,1) + abs(par.bc(par.k,extXs(:,av)) - evalFieldQBF(par,extXs(:,av),c1,0) )/vSecKind.avm;
    vSecKind.errExt(1,2) = vSecKind.errExt(1,2) + abs(par.bc(par.k,extXs(:,av)) - evalFieldQBF(par,extXs(:,av),c2,1) )/vSecKind.avm;
end
vSecKind.errExt = vSecKind.errExt/mean(abs(asdf(1,:)));
figure; plot(taus, [abs(asdf(2,:)); abs(asdf(2,:) -asdf(1,:)); abs(asdf(3,:) -asdf(1,:))] );
vSecKind  % The error on the boundary conditions has not been adjusted to Neumann BC
% v = validate(A1,A2, par, struct('field', zeros(70))); plotVal(v)
return
%% test
for ix = 1:100:par.N
    figure; plot(par.colltau, [real(transpose(A1(ix,:)).*c1), imag(transpose(A1(ix,:)).*c1)]); title(num2str(ix))
end
figure; plot(par.colltau, [real(c1), imag(c1)]); title('Solution vector c1');


%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

par = getObst(1);
Tcor = 0.02;
percDecay = 1;
% ks = [2^7, fzero(@(x) besselj(0,x), 130)];
% sk = 2*besselroots(100,1);

avm = 100;
% vAac = struct('perc', zeros(1,2), 'errSol', zeros(1,2), 'conds', zeros(1,2), ...
%     'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(1,2) );

%% Calculations (close to and) at resonance
% for ki = 1:kl
printtoc = 5;
% par.k = ks(ki);
par.k = 2*5.830649350672455e+01; % Two times the first zero of J_{50}(x)
par.N = ceil(par.ppw*par.k);
par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
par.colltau = par.t(1:par.N);

%% Computations for the full matrix
A1 = zeros(par.N);
tic;
prevToc = toc;
for i = 1:par.N
    % 	parfor i=1:par.N % Possibly use this instead of a sequential loop
    A1(i,:) = collRowQBF(i,par);
end

%% Compute the correlations
% if ki == 1
b = par.bc(par.k,par.par(par.colltau));
c1 = A1\b;
[R, sigma] = calcCorr(par, c1, Tcor, percDecay, [], A1);
colLow = par.colltau;
% end

%% Computations for standard asymptotically compressed matrix
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

%% Validation

vAac = validate(A1,A2,par,struct('perc', zeros(1,2), 'errSol', zeros(1,2), 'conds', zeros(1,2), ...
    'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(1,2) ),1)
% vAac = validate(A1,A2,par,vAac,ki)

% end

return

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


%% Cubic basis functions
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

return


%% Test resonance
clearvars

par = getObst(1);
printtoc = 5;
% J_0(ka) = 0
% par.k = fzero(@(x) bessely(0,x), 128)*2;
% par.k = fzero(@(x) besselj(0,x), 1.1)*2;
% par.k = fzero(@(x) besselj(1,x), 2.1)*2;
% par.k = fzero(@(x) besselj(1,x), 2.1)*2.001;
% par.k = fzero(@(x) besselj(0,x), 1.1)*2.001;
% par.k = fzero(@(x) bessely(0,x), 1.1)*(1.1);%+sqrt(eps));
moBesl = 22;
bnd = besselroots(moBesl,1)*2;
albesl = zeros(moBesl);
for mob = 0:moBesl-1
    bnd = [bnd; besselroots(mob,moBesl)*2];
    albesl(:,mob+1) = besselroots(mob,moBesl)*2;
end
ks = sort(bnd(bnd <= bnd(1)) );


coA = nan*ks;
nmis = nan*ks;
izxs = nan*ks;
for ki = 1:length(ks)
    par.k = ks(ki);
    par.N = ceil(par.ppw*par.k);
    par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
    par.colltau = par.t(1:par.N);
    
    A1 = zeros(par.N);
    tic;
    prevToc = toc;
    for i = 1:par.N
        if (toc-prevToc > printtoc)
            prevToc = toc;
            display(['k = ', num2str(par.k), ', ' num2str( (i-1)/par.N,'%7.3f'), ' = part of matrix, estim sec left = ' num2str(toc*(par.N-i+1)/(i-1)) ])
        end
        % 	parfor i=1:par.N % Possibly use this instead of a sequential loop
        A1(i,:) = collRowQBF(i,par);
    end
    coA(ki) = cond(A1);
    [nmi, izx] = find(ks(ki) == albesl);
    nmis(ki) = nmi;
    izxs(ki) = izx;
    display([num2str(ki) '=ki, kl=', num2str(length(ks)), ', cond=' num2str(coA(ki)), ', k=' num2str(ks(ki)), ...
        ', besOrder=', num2str(nmi), ', ixZero = ' num2str(izx)])
end
