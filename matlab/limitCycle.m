% Compute the solution of multiple scattering obstacles using a known phase.

%% Initialising
clearvars
% close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

ks = 2^9;
obsts = 5;
% obsts = 11;
% rx = 0.1;
% rx = 0.3;
% rx = 0.7;
% rx = 1.4;
% rx = 4.4;
bc = 1;
% bcsh = pi/2; % shift in boundary condition

printtoc = 3;
kl = length(ks);
nbOb = length(obsts);

avm = 100; % Number of random taus to average BC over
v = struct('avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(nbOb*kl,2), 'errSol', zeros(nbOb*kl,2),  ...
    'errInt', zeros(nbOb*kl,2), 'timeA', zeros(nbOb*kl,4), 'ks', ks); %, 'field', zeros(70) );

%% Computations
for oi = 1:length(obsts)
    obstacle = obsts(oi);
    start = now;
    
    for ki = 1:kl
        idx = (oi-1)*kl+ki;
        if obstacle == 11
            par = getObst(obstacle, rx); % Reset par
        else
            par = getObst(obstacle); % Reset par
            rx = 0.5;
        end
        par.k = ks(ki);
        par.N = 0;
        par.r = zeros(2,length(par.obsts)); % ranges
        for obst = 1:length(par.obsts)
%             par.obsts(obst).ppw = par.obsts(obst).ppw*2; % Test: increase sampling rate
            par.obsts(obst).k = par.k;
            par.obsts(obst).N = par.obsts(obst).ppw*par.k;
            par.r(1,obst) = par.N+1;
            par.N = par.N + par.obsts(obst).N;
            par.r(2,obst) = par.N;
            
            par.obsts(obst).t = linspace(0,1,par.obsts(obst).N+1); % The knots of the periodic spline;
            par.obsts(obst).colltau = par.obsts(obst).t(1:par.obsts(obst).N);
            par.obsts(obst).hs = (par.obsts(obst).colltau(2) -par.obsts(obst).colltau(1) )/2;
        end
        
        if bc == 4
            par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/2-bcsh);sin(-pi/2-bcsh)]);
        elseif bc == 3
            par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/2);sin(-pi/2)]);
        elseif bc == 2
            par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/4);sin(-pi/4)]);
        elseif bc == 1
            par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(0);sin(0)]);
        end
        
        overs1 = 1;
        
        %% Computating full solution
        A1 = zeros(overs1*par.N, par.N); 
        b1 = zeros(overs1*par.N,1);
        tic
        prevToc = toc;
        for i = 1:par.N
            if (toc-prevToc > printtoc)
                prevToc = toc;
                display([num2str(par.k), '=k, ' num2str( (i-1)/par.N,'%7.3f'), '=A1%, est. # sec. left for A1=' num2str(toc*(par.N-i+1)/(i-1)) ])
            end
            obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
            collx = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) );
            
            collxos = zeros(2,overs1);
            for os = 1:overs1
                collxos(:,os) = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) + par.obsts(obsin).hs*2/overs1*(os-1));
            end
            for obst = 1:length(par.obsts)
                if obst == obsin
                     for os = 1:overs1
                        A1(overs1*i-overs1+os,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N, [], [], ...
                            par.obsts(obst).hs*2/overs1*(os-1) );
                        b1(overs1*i -overs1+os) = par.bc(par.k, collxos(:,os));
                     end
                else
                    for os = 1:overs1
                        A1(overs1*i-overs1+os,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxos(:,os) );
                    end
                end
            end
        end
        v.timeA(idx,1) = toc;
        c1 = A1\b1; 
        % [R, sigma,obbounds] = calcCorr(par, c1, Tcor, percDecay); 
        % [R, sigma,obbounds] = calcCorr(par, c1, 0.1, 1, [10, now], A1); 
                
        %% Some plots
%         figure; spectrogram(c1, round(length(c1)/16),[],[], 1/(par.obsts(1).colltau(2)-par.obsts(1).colltau(1)), 'centered'); title(['Spectrogram ' ])%num2str(bcsh)])
%         figure; plot( [real(c1) imag(c1)]); legend('Re(c1)','Im(c1)')
%         v = validate(A1,nan*A1,par,v,idx)

    end
    display([num2str(oi) ' = oi, ki = ' num2str(ki) ', now is ' datestr(now) ', expected end ' datestr(start + ...
        (now-start)*sum(ks.^2)/sum(ks(1:ki).^2)*( length(obsts)-oi + 1) )  ]);
end


%% Ray tracing by partitioned block 
A11 = A1(1:size(A1,1)/2, 1:size(A1,2)/2);
A21 = A1(1:size(A1,1)/2, (1+size(A1,2)/2):end); % Action of density of obst 2 on location of obst 1
A12 = A1((size(A1,1)/2+1):end, 1:size(A1,2)/2);
A22 = A1((size(A1,1)/2+1):end, (1+size(A1,2)/2):end);

c1true1 = c1(1:length(c1)/2);
c1true2 = c1(length(c1)/2+1:end);

nbrefl = 1; %8;
bsl1 = zeros(size(A1,1)/2, nbrefl);
bsl2 = zeros(size(A1,1)/2, nbrefl);
csl1 = zeros(size(A1,2)/2, nbrefl);
csl2 = zeros(size(A1,2)/2, nbrefl);

bsl1(:,1) = b1(1:size(A1,1)/2);
bsl2(:,1) = b1((size(A1,1)/2+1):end);
csl1(:,1) = A11\bsl1(:,1);
csl2(:,1) = A22\bsl2(:,1);

resid = cell(2,nbrefl);
resid{1,1} = c1true1 - sum(csl1,2);
% figure;
% plot(par.obsts(1).colltau, real(csl1(:,1)) ); title('csl1 refl')
% hold on

for refl = 2:nbrefl
    bsl1(:,refl) = -A21*csl2(:,refl-1);
    bsl2(:,refl) = -A12*csl1(:,refl-1);
    csl1(:,refl) = A11\bsl1(:,refl);
    csl2(:,refl) = A22\bsl2(:,refl);
    
    resid{1,refl} = c1true1 - sum(csl1,2);
%     plot(par.obsts(1).colltau, real(csl1(:,refl)));
end
% legend(num2str((1:nbrefl)'));


%% Eigenvalues big reflection matrix

bigm = A11\A21*(A22\A12);
display(['Starting EVD of matrix of size ' num2str(size(bigm))]);
[Vb, Db, Wb] = eig(bigm);
display(['Ended EVD of matrix of size ' num2str(size(bigm))]);
% [V1, D1, W1] = eig(A1);
% figure; semilogy(abs(diag(Db)))
figure; plot(par.obsts(1).colltau, [real(Vb(:,1)) imag(Vb(:,1))]); legend('Re Vb', 'Im Vb');


c1nor = sqrt(sum(abs(csl1).^2, 1));
c2nor = sqrt(sum(abs(csl1).^2, 1));
ratcnor = c1nor(2:end)./c1nor(1:end-1).*c2nor(2:end)./c2nor(1:end-1)
% Modulus also for negative angles
[real(Db(1,1)) imag(Db(1,1)) abs(Db(1,1)) mod(angle(Db(1,1)),2*pi) mod(2*1*par.k,2*pi)]



%% Calculated tau2 (D261) for a circle

c1x = -0.5; c1y = -0.5; c2x = -0.5; c2y = 1.5; r1 = rx; r2 = 0.5;

% tau2 = (-0.5-1.5+0.5*sin(2*pi*par.obsts(1).colltau) )./sqrt((-0.5-1.5)^2+2*0.5*sin(2*pi*par.obsts(1).colltau)*(-0.5-1.5) + 0.5^2);
tau2 = (c1y-c2y+r2*sin(2*pi*par.obsts(1).colltau) )./sqrt((c1y-c2y)^2+2*r1*sin(2*pi*par.obsts(1).colltau)*(c1y-c2y) + r1^2);
% tau2 = asin(tau2)/2/pi+1;
% tau2 = (-1).^(par.obsts(1).colltau > 0.25).*asin(tau2)/2/pi+1 -0.5*(par.obsts(1).colltau > 0.25);
tau2 = (-1).^(abs(par.obsts(1).colltau-0.5) < 0.25).*asin(tau2)/2/pi+1 -0.5*(abs(par.obsts(1).colltau-0.5) < 0.25);
% figure; plot(par.obsts(1).colltau, tau2); title('tau2 D261')
% tau2 = 0.75 + (tau2-0.75)/2;
figure; plot(par.obsts(1).colltau, tau2);
0.25+max(tau2)
asin(0.5/2)/2/pi

phider2 = -( (-0.5-(-0.5)+rx*(cos(2*pi*par.obsts(1).colltau) - cos(2*pi*tau2) ) ).*rx.*2.*pi.*sin(2*pi*tau2) -0.5*2*pi*cos(2*pi*tau2).*...
    (-0.5-1.5 + 0.5*sin(2*pi*par.obsts(1).colltau) - 0.5*sin(2*pi*tau2) ) )./sqrt( (-0.5-(-0.5)+rx*(cos(2*pi*par.obsts(1).colltau) - cos(2*pi*tau2) ) ).^2 +...
    (-0.5-1.5 + 0.5*sin(2*pi*par.obsts(1).colltau) - 0.5*sin(2*pi*tau2) ).^2 );
% figure; plot(par.obsts(1).colltau, phider2);
% figure; plot(tau2, phider2);
% figure; plot(par.obsts(1).colltau, cumsum(phider2));
phi = sqrt( (c1x-c2x+r1*cos(2*pi*par.obsts(1).colltau)-r1*cos(2*pi*tau2) ).^2 + (c1y-c2y+r2*sin(2*pi*par.obsts(1).colltau)-r2*sin(2*pi*tau2) ).^2 );
% figure; plot(par.obsts(1).colltau, phi); title('phase');


%% Extract phase from angle
collsignal = par.obsts(1).colltau;
% signal = Vb(:,1);
signal = Vb1k9;
% signal = Vb(:,1)./abs(Vb(:,1));
% figure; plot(collsignal, [real(signal) imag(signal)]);
l = length(signal)/2;
ph = zeros(size(signal));
tryph = transpose(sqrt( sum( (par.obsts(obst).par(par.obsts(obst).colltau) - repmat(par.obsts(3-obst).par((5-2*obst)/4), 1, ...
    length(par.obsts(obst).colltau) ) ).^2, 1) ) );
% tryph = transpose(2.1618/2*(cos(4*pi*par.obsts(obst).colltau)-1));
% tryph = arrayfun( @(tau1) tau1, transpose(par.obsts(obst).colltau));
for i = l:-1:1
    ph(i) = angle(signal(i));
    if (i ~= l) && (abs(ph(i) - ph(i+1)) > 5)
        ph(i) = ph(i) -2*pi*round( (ph(i)-ph(i+1))/2/pi);
    end
end
for i = l:length(signal)
    ph(i) = angle(signal(i));
    if (i ~= l) && (abs(ph(i) - ph(i-1)) > 5)
        ph(i) = ph(i) -2*pi*round( (ph(i)-ph(i-1))/2/pi);
    end
end
closest = find(abs(collsignal-0.25) == min(abs(collsignal-0.25)));
ph = ph + 1*ks(ki) - ph(closest);
% figure; plot(collsignal, [ph/ks(ki) tryph]); legend('phase from angle', 'try phase');
figure; plot(collsignal, [ph/ks(ki) tryph transpose(phi)]); legend('phase from angle', 'try phase', '\phi from closest \tau_2');
figure; plot(collsignal, ph/ks(ki)-(tryph+transpose(phi))/2); legend('phase from angle - average');
norm(ph(500:2000)/ks(ki)-(tryph(500:2000)+transpose(phi(500:2000)))/2)/norm(ph(500:2000)/ks(ki)) % for ks = 2^9
% figure; plot(collsignal(1:end-1), (ph(2:end)-ph(1:end-1))/(collsignal(2)-collsignal(1))/ks(ki)); legend('der phase from angle');
% figure; plot(collsignal, [real(signal.*exp(1i*par.k*tryph)) imag(signal.*exp(1i*par.k*tryph))]); legend('real LF

% return
%% Determine tau2 from ph as a distance
t2d = zeros(l,3);
disd = @(t1, t2) sqrt(r1^2*(cos(2*pi*t1)-cos(2*pi*t2)).^2 + (c1y-c2y+r1*sin(2*pi*t1)-r1*sin(2*pi*t2)).^2);
for i=1:l
    t2d(i,1) = fsolve(@(t2) ph(i)/ks - disd(collsignal(i), t2), 0.75, optimoptions('fsolve', 'Display', 'none'));
    t2d(i,2) = ph(i)/ks - disd(collsignal(i), t2d(i,1));
%     t2d(i,3) = t2d(i,2) - phi(tau2-1/2);
end

for i = 1:l
    [~, p] = min(abs(t2d(i,1) - 1/2 - collsignal));
    t2d(i,3) = (phi(p+1) - phi(p))/(collsignal(2) - collsignal(1)) +(r1^2*(cos(2*pi*collsignal(i)) - cos(2*pi*t2d(i)))*2*pi*sin(2*pi*t2d(i)) +...
        (c1y-c2y +r2*sin(2*pi*collsignal(i)) - r2*sin(2*pi*t2d(i,1)))*(-r1*2*pi)*cos(2*pi*t2d(i)) )/disd(collsignal(i), t2d(i,1));
end
asdftau2 = (c1y-c2y+r2*sin(2*pi*collsignal(1:l)) )./sqrt((c1y-c2y)^2+2*r1*sin(2*pi*collsignal(1:l) )*(c1y-c2y) + r1^2);
figure; plot(collsignal(1:l), t2d(:,1)); title('tau2 via ph as dist')
figure; plot(collsignal(1:l), [sin(2*pi*t2d(:,1)), asdftau2']);
asdftau2 = (-1).^(abs(par.obsts(1).colltau(1:l)-0.5) < 0.25).*asin(asdftau2)/2/pi+1 -0.5*(abs(par.obsts(1).colltau(1:l)-0.5) < 0.25);
figure; plot(collsignal(1:l), [t2d(:,1), asdftau2']);

pm0 = [(c1y-c2y), r1, ((c1y-c2y)^2+r1^2), 2*r1*(c1y-c2y)]; %[rand(4,1)];
[pm, fv] = fsolve( @(pm) sin(2*pi*t2d(:,1))-(pm(1)+pm(2)*sin(2*pi*collsignal(1:l)'))./sqrt(pm(3) + pm(4)*sin(2*pi*collsignal(1:l)')), pm0);

pm0 = [(c1y-c2y), r1, ((c1y-c2y)^2+r1^2), 2*r1*(c1y-c2y), 0, 0];
[pm, fv] = fsolve( @(pm) sin(2*pi*t2d(:,1))-(pm(1)+pm(2)*sin(2*pi*collsignal(1:l)')+pm(5)*sin(2*pi*collsignal(1:l)').^2)./...
    sqrt(pm(3) + pm(4)*sin(2*pi*collsignal(1:l)')+pm(6)*sin(2*pi*collsignal(1:l)').^2), pm0);

%% 'exact' integral
th1f = @(t) asin( ( -(2*r1*(c1y-c2y)-2*r1*(c1y-c2y)*sin(2*pi*t).^2) + sqrt( (2*r1*(c1y-c2y)-2*r1*(c1y-c2y)*sin(2*pi*t).^2).^2 ...
    -4*r1^2*( (c1y-c2y)^2*(1-sin(2*pi*t).^2)-r1^2*(sin(2*pi*t)).^2)) )/2/r1^2);
% figure; plot(collsignal, [real(th1f(collsignal)/2/pi)' imag(th1f(collsignal)/2/pi)']);
tau1f = @(t) (-1).^(abs(t-0.5) < 0.25).*th1f(t)/2/pi+1 +0.5*(abs(t-0.5) < 0.25)-1;
figure; plot(collsignal, [real(tau1f(collsignal))' imag(tau1f(collsignal))']);

dist = @(t) sqrt(r1^2*(cos(2*pi*tau1f(t))-cos(2*pi*t)).^2 + (c1y-c2y+r1*sin(2*pi*tau1f(t))-r1*sin(2*pi*t)).^2);
integr = @(t) (r1*(cos(2*pi*tau1f(t))-cos(2*pi*t)).*2.*pi.*sin(2*pi*t) +(c1y-c2y+r1*sin(2*pi*tau1f(t))-r1*sin(2*pi*t)).*(-2*r1*pi).*cos(2*pi*t) )./dist(t);
t2s = linspace(min(tau2), max(tau2), 20);
phs = 0*t2s;
vlg = 0*t2s;
for i = 1:length(t2s)
%     phs(i) = quadgk( @(t) integr(t), 0.75, t2s(i) ) + 2*1;
    phs(i) = 1 - quadgk( @(t) integr(t), 0.75, t2s(i) )*sqrt(2);
end
for i = 1:length(t2s)
    [~, m] = min(abs(t2s-th1f(t2s(i))));
    vlg(i) = dist(t2s(i)) + phs(i) - phs(m);
end
% figure; plot(t2s, phs); title('phase via integration over tau2');
figure; plot(collsignal, [ph/ks(ki) tryph transpose(phi)]); 
hold on; plot(t2s-0.5, phs); legend('phase from angle', 'try phase', '\phi from closest \tau_2', 'int over tau2');



%% Project eigenvectors of bigm onto assumed phase
obst = 1;
% phbigm = exp(1i*par.k*sqrt( sum( (par.obsts(obst).par(par.obsts(obst).colltau) - repmat(par.obsts(3-obst).par((2*obst-1)/4), 1, ... % Wrong
%     length(par.obsts(obst).colltau) ) ).^2, 1) ) ); 
% phbigm = exp(1i*par.k*sqrt( sum( (par.obsts(obst).par(par.obsts(obst).colltau) - repmat(par.obsts(3-obst).par((5-2*obst)/4), 1, ...
%     length(par.obsts(obst).colltau) ) ).^2, 1) ) );
% % phbigm = exp(1i*par.k*sqrt( sum( (par.obsts(obst).par(par.obsts(obst).colltau) - repmat(par.obsts(3-obst).par(arcsin(sqrt(()/()))*), 1, ...
% %     length(par.obsts(obst).colltau) ) ).^2, 1) ) );
phbigm = exp(1i*par.k*phi);
% phbigm = exp(1i*par.k*(phi+sqrt( sum( (par.obsts(obst).par(par.obsts(obst).colltau) - repmat(par.obsts(3-obst).par((5-2*obst)/4), 1, ...
%     length(par.obsts(obst).colltau) ) ).^2, 1) ) )/2 );
lfbm = Vb(:,1)./transpose(phbigm);
figure; 
plot(par.obsts(1).colltau, [real(lfbm) imag(lfbm)])
title(['LF behaviour if first eigenvector for rx = ' num2str(rx)]);
% figure; plot(par.obsts(1).colltau, abs(lfbm)); title('abs lfbm')
% figure; plot(par.obsts(1).colltau, angle(lfbm)); title('angle lfbm')
% plot cos(2*pi*tau) /(arctanh((sin(2*pi*tau)-1)/sqrt(-cos(2*pi*tau)^2)) - arctanh((sin(2*pi*tau)+1)/sqrt(-cos(2*pi*tau)^2)) )for tau in [0.25,0.75]


%% Nonlinear solve
% nr = 100;
% nr = 1e5;
% nr = 1e4;
nr = 1e3;
cl1 = linspace(0,1/2,nr);
c1x = -0.5; c1y = -0.5; c2x = -0.5; c2y = 1.5;
% c1x = -0.5; c1y = -0.5; c2x = -0.5; c2y = 11.5;
r1 = rx; r2 = 0.5;
h = cl1(2)-cl1(1);
% sel = @(y,z) y(z);
% d = -2*norm(par.obsts(1).par(0.25) - par.obsts(2).par(0.75));
d = 0;
% clos = @(phi, tau2, shft) arrayfun(@(t2) phi(find( abs(cl1-t2) == min(abs(cl1-t2)))+shft), tau2);
clos = @(phi,tau2, shft) interp1(cl1-shft*h, phi, tau2);
dist = @(tau2) sqrt( (c1x-c2x+r1*cos(2*pi*cl1)-r1*cos(2*pi*tau2) ).^2 + (c1y-c2y+r2*sin(2*pi*cl1)-r2*sin(2*pi*tau2) ).^2 );

F = @(x) [ (x(1:nr)-dist(x(nr+1:end))-clos(x(1:nr), x(nr+1:end)-1/2, 0)-d ), ...
    ((clos(x(1:nr),x(nr+1:end)-1/2, 1)-clos(x(1:nr),x(nr+1:end)-1/2, 0))/h + ...
    2*pi*(r1*sin(2*pi*x(nr+1:end) ).*(c1x-c2x+r1*cos(2*pi*cl1)-r1*cos(2*pi*x(nr+1:end))) ...
    -r2*cos(2*pi*x(nr+1:end)).*(c1y-c2y+r2*sin(2*pi*cl1)+r2*sin(2*pi*x(nr+1:end))) )./dist(x(nr+1:end) ) ) ];

x0 = [sqrt( sum( (par.obsts(1).par(cl1) - repmat(par.obsts(2).par(0.75), 1, length(cl1) ) ).^2, 1) ), 0.75*ones(size(cl1))];
% x0 = [rand(1,nr), 0.6+0.2*rand(1,nr)];
opto = optimoptions('fsolve', 'Display', 'iter-detailed', 'FunctionTolerance', 1e-2, 'OptimalityTolerance', 1e-3, 'StepTolerance', 1e-12);
% opto = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter-detailed', 'FunctionTolerance', 1e-2*norm(F(x0)),...
%     'OptimalityTolerance', 1e-3, 'StepTolerance', 1e-3);
[X, FVAL] = fsolve(F, x0, opto);

trph = sqrt( sum( (par.obsts(obst).par(cl1) - repmat(par.obsts(2).par(0.75), 1, nr) ).^2, 1) );
figure; plot(cl1, X(1:nr)- trph); title('diff \phi')
figure; plot(cl1, X((nr+1):end)); title('\tau_2')
% figure; plot(cl1, [X(1:100); trph]); title('\phi')

