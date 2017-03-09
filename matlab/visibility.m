% Use the visibility criterion to determine which regions to include in the windows for the near-inclusion obstacle.
%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

par = getObst(3);

printtoc = 5;
par.k = 2.^9;
par.N = par.ppw*par.k;
par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
par.colltau = par.t(1:par.N);

%% Plot of the obstacle and visible region in the parameter domain
mnr = 200;
ts = linspace(0,1,mnr)';
parametr = par.par(ts');
fs = 25;
%$[\min(0.365-T,\tau-0.02-T), \max(0.665+T,\tau+0.02+T)]$ for $\tau \in [0.3,0.72]$, where the length of the decaying part $T=0.1$
leftb = 0.365;
rightb = 0.665;
T = 0.1;

figure;
plot(parametr(1,:),parametr(2,:),'k', 'LineWidth',3);
hold on;
tsSight = ts((ts >= leftb) & (ts <= rightb) )';
pSight = par.par(tsSight)-0.008*par.normal(tsSight);
plot(pSight(1,:),pSight(2,:),'b','LineWidth',5);

tsRoloff = ts((ts >= leftb -T) & (ts < leftb) )';
pRol = par.par(tsRoloff)-0.008*par.normal(tsRoloff);
plot(pRol(1,:),pRol(2,:),'b:','LineWidth',5);

tsRoloff = ts((ts > rightb) & (ts <= rightb +T))';
pRol = par.par(tsRoloff)-0.008*par.normal(tsRoloff);
plot(pRol(1,:),pRol(2,:),'b:','LineWidth',5);

col = par.par(0.5);
up = (par.par(leftb)-col)*1.56+col;
down = (par.par(rightb)-col)*1.5+col;
plot([up(1) col(1) down(1)], [up(2) col(2) down(2)], 'k--');
te = linspace(0,0.95,20);
spar = par.par(te);
plot(spar(1,:), spar(2,:), 'r*','LineWidth',5);
labels = cellstr( num2str(te.') );
h = text(spar(1,:), spar(2,:), labels, 'VerticalAlignment','bottom', ...
	'HorizontalAlignment','right');
set(h,'FontSize',fs);
set(gca,'FontSize',fs);
hold off;

%% Computations for the full matrix
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


for doCompr = 0:1
    %% Computations for compressed matrix
    A2 = zeros(par.N);
    tic;
    prevToc = toc;
    for i = 1:par.N
        if (toc-prevToc > printtoc)
            prevToc = toc;
            display([num2str(T,'%8.4f') '=T, ' num2str(i/par.N,'%7.3f') 'A2%, est. sec left=' num2str(toc*(par.N-i)/i)])
        end
        tc = par.colltau(i);
        if (tc < 0.2)
            bounds = [];
        elseif tc < 0.3
            if doCompr
                bounds = [tc-0.02-T; tc+0.02+T; T; T];
            else
                bounds = [];% Use this to do no compression for these collocation points either
            end
        elseif tc < 0.72
            bounds = [min(leftb-T,tc-0.02-T); max(rightb+T,tc+0.02+T); T; T];
            %$[\min(0.365-T,\tau-0.02-T), \max(0.665+T,\tau+0.02+T)]$ for $\tau \in [0.3,0.72]$, where the length of the decaying part $T=0.1$
        elseif tc < 0.8
            if doCompr
                bounds = [tc-0.02-T; tc+0.02+T; T; T];
            else
                bounds = [];
            end
        else
            bounds = [];
        end
        A2(i,:) = windRow(i,par,bounds);
    end % loop over row-indices i
    
    %% Show the results
    plotVal(struct(),A2) % Show the structure of the compressed matrix
    avm = 100;
    display([num2str(par.k) ' = k, doCompr = ' num2str(doCompr)]);
    v = validate(A1,A2,par,struct('perc', zeros(1,2), 'errSol', zeros(1,2), 'conds', zeros(1,2), ...
        'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(1,2) ),1)
end



%% Two circles: initialisation of the calculation
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

par = getObst(5);
Tcor = 0.15;
percDecay = 1;
printtoc = 0.5;
par.k = 2.^8;

mti = 0;
avm = 100; % Number of random taus to average BC over
v = struct('mti', mti, 'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(1,2+mti),...
    'nnz', zeros(1,2), 'perc', zeros(1,2), 'errSol', zeros(1,2+mti), ...
    'errBCcol', zeros(1,2+mti), 'compresErr', zeros(1,2),  'errInt', zeros(1,2+mti), ...
    'timeSol', zeros(1,2+mti), 'nbIter', zeros(1,mti), 'timeA', zeros(1,4), 'ks', par.k);

start = now;
startk = now;
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

%% Computating the full solution
A1 = zeros(par.N); tic; prevToc = toc;
for i = 1:par.N
    if (toc-prevToc > printtoc)
        prevToc = toc;
        display([num2str(i/par.N,'%7.3f') '=A1%, now=' datestr(now) ', est. # sec. left for A1=' num2str(toc*(par.N-i)/i) ]);
    end
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
v.timeA(1,1) = toc;

% The solution is needed for computing the correlations
b = zeros(par.N,1);
for obst = 1:length(par.obsts)
    b(par.r(1,obst):par.r(2,obst)) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).colltau));
end
c1 = A1\b;

if 0
    [R, sigma] = calcCorr(par, c1, Tcor, percDecay, [5, now()], A1); % The approximation of R used to make figure
    figure;
    pcolor(log(max(abs(R),1e-2)) );
    hold on; shading interp;
    xlabel('n'); ylabel('m');
    hh = colorbar();
    ylabel(hh,'$\log(\max(|R_{m,n}|,10^{-2}))$', 'interpreter','latex', 'FontSize',20);
    set(gca,'YDir','rev')
end

%% Compressed matrix
% comprMeth = 1; % then only compression in the illuminated region by the source, 
% comprMeth = 2; % then correct compression in generalised illuminated region by both source and other circle,
comprMeth = 3; % then correct illuminated and compression in coupling matrices

A2 = zeros(par.N);
tic;
prevToc = toc;
for i = 1:par.N
    if (toc-prevToc > printtoc)
        prevToc = toc;
        display([num2str(i/par.N,'%7.3f') '=A2%, now=' datestr(now) ', est. # sec. left for A2=' ...
            num2str(toc*(par.N-i)/i) ', tmp comprErr=' num2str(norm(A2(1:(i-1),:)*c1-b(1:(i-1)))/norm(b(1:(i-1))  )) ]);
    end
    obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
    tc = par.obsts(obsin).colltau(i-par.r(1,obsin)+1);
    collx = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) );
    
    for obst = 1:length(par.obsts)
        if (obst == obsin) 
            bounds = [];
            pt = par.obsts(obsin).par(tc);
            if ((comprMeth == 1) && (pt(1) < -0.5-Tcor)) || ( (comprMeth >= 2) &&...
                   ( ((obsin == 1) && (tc > 0.27) && (tc < 0.48)) || ((obsin == 2) && (tc > 0.52) && (tc < 0.73)) ) )
                bounds = [tc-2*Tcor; tc+2*Tcor; percDecay*Tcor; percDecay*Tcor];
            end
            A2(i,par.r(1,obsin):par.r(2,obsin)) = windRow(i-par.r(1,obsin)+1, par.obsts(obsin), bounds);
        elseif obsin == 2
            if tc < 0.5
                bounds = [0; 0.5];
            elseif tc < 3/4
                g2g2 = acos(4*sin(2*pi*tc)*cos(2*pi*tc)+cos(2*pi*tc) -1/2*sqrt(cos(2*pi*tc)^2*(8*sin(...
                    2*pi*tc)+2)^2 -4*(15*sin(2*pi*tc)^2+8*sin(2*pi*tc)+1)) );
                if imag(g2g2) ~= 0
                    g2g2 = -inf;
                end
                bounds = [max(g2g2, asin((8+2*sin(2*pi*tc) + sqrt(64+32*sin(2*pi*tc)+4*sin(2*pi*tc)^2 ...
                    -4*sin(2*pi*tc)^2*(17+8*sin(2*pi*tc))))/(34+16*sin(2*pi*tc))) ); ...
                    (pi-asin((8+2*sin(2*pi*tc) -sqrt(64+32*sin(2*pi*tc)+4*sin(2*pi*tc)^2 ...
                    -4*sin(2*pi*tc)^2*(17+8*sin(2*pi*tc))))/(34+16*sin(2*pi*tc))) )]/2/pi;
            else
                g2g2 = acos(4*sin(2*pi*tc)*cos(2*pi*tc)+cos(2*pi*tc) +1/2*sqrt(cos(2*pi*tc)^2*(8*sin(...
                    2*pi*tc)+2)^2 -4*(15*sin(2*pi*tc)^2+8*sin(2*pi*tc)+1)) );
                if imag(g2g2) ~= 0
                    g2g2 = inf;
                end
                bounds = [asin((8+2*sin(2*pi*tc) -sqrt(64+32*sin(2*pi*tc)+4*sin(2*pi*tc)^2 ...
                    -4*sin(2*pi*tc)^2*(17+8*sin(2*pi*tc))))/(34+16*sin(2*pi*tc))); ...
                    min(g2g2, pi-asin((8+2*sin(2*pi*tc) + sqrt(64+32*sin(2*pi*tc)+4*sin(2*pi*tc)^2 ...
                    -4*sin(2*pi*tc)^2*(17+8*sin(2*pi*tc))))/(34+16*sin(2*pi*tc))))]/2/pi;
            end
            if comprMeth >= 3
                bounds = [0; 0.5]; % Ignore the bounds from the part of obst 1 that can be seen from obst 2
                decay = repmat(Tcor*percDecay,1,size(bounds,2));
                bounds = [(bounds + [-decay; decay]/percDecay); [1;1]*decay];
                A2(i,par.r(1,obst):par.r(2,obst)) = windRow('error', par.obsts(obst), bounds,[],[],collx);
            else
                A2(i,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N,  [], collx);
            end
        else % obsin = 1, obst = 2
            if comprMeth >= 3
                bounds = [0.5; 1];
                decay = repmat(Tcor*percDecay,1,size(bounds,2));
                bounds = [(bounds + [-decay; decay]/percDecay); [1;1]*decay];
                A2(i,par.r(1,obst):par.r(2,obst)) = windRow('error', par.obsts(obst), bounds,[],[],collx);
            else
                A2(i,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N,  [], collx);
            end
        end
    end
end
v.timeA(1,2) = toc;
v = validate(A1,A2,par,v,1)

figure; spy(A2);
set(gca,'XTick',0:1e3:5e3);
set(gca,'tickdir','out'); 
ylabel('m'); 
xlabel('l'); 
set(gca,'Fontsize', 22)

