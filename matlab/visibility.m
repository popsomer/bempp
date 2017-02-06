% Use the visibility criterion to determine which regions to include in the windows for the near-inclusion obstacle.
%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

par = getObst(3);

printtoc = 5;
par.k = 2.^8;
par.N = par.ppw*par.k;
par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
par.colltau = par.t(1:par.N);

%% Plot of the obstacle and visible region in the parameter domain
mnr = 200;
ts = linspace(0,1,mnr)';
parametr = par.par(ts');
fs = 25;


leftb = 0.27;
rightb = 0.73;

% Do not take T large here when bounds(1:2,:) are fixed
T = 0.02;

figure;
plot(parametr(1,:),parametr(2,:),'k', 'LineWidth',3);
hold on;
tsSight = ts((ts >= leftb-T) & (ts <= rightb) )';
pSight = par.par(tsSight)-0.008*par.normal(tsSight);
plot(pSight(1,:),pSight(2,:),'b','LineWidth',5);

tsRoloff = ts((ts >= lefte) & (ts < leftb +T) )';
pRol = par.par(tsRoloff)-0.008*par.normal(tsRoloff);
plot(pRol(1,:),pRol(2,:),'b:','LineWidth',5);
tsRoloff = ts((ts > rightb -T) & (ts <= righte))';
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
		bounds = [min(tc-T*1.1,0.27);max(0.25,tc+T*1.1);T;T];
	elseif tc < 0.7
		bounds = [leftb-T; rightb+T; T; T];
	elseif tc < 0.8
		bounds = [min(0.75,tc-T*1.1); max(tc+T*1.1,0.73);T;T];
	else
		bounds = [];
	end
	A2(i,:) = windRow(i,par,bounds);
end % loop over row-indices i

%% Show the results
plotVal(struct(),A2) % Show the structure of the compressed matrix
% v = validate(A1,A2,par,struct('perc', zeros(1,2), 'errSol', zeros(1,2), 'conds', zeros(1,2) ),1)
avm = 100;
v = validate(A1,A2,par,struct('perc', zeros(1,2), 'errSol', zeros(1,2), 'conds', zeros(1,2), ...
    'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(1,2) ),1)
