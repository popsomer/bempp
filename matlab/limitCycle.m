% Two circlesS

% %Compute the solution of multiple scattering obstacles using a known phase.

%% Initialising
clearvars
% close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

ks = 2^9;
% ks = 2^7;
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
% load Vb1k9; 
signal = Vb1k9;
% signal = Vb(:,1)./abs(Vb(:,1));
% figure; plot(collsignal, [real(signal) imag(signal)]);
l = length(signal)/2;
ph = zeros(size(signal));
obst = 1;
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
% figure; plot(collsignal, [ph/ks(ki) tryph transpose(phi)]); legend('phase from angle', 'try phase', '\phi from closest \tau_2');
% figure; plot(collsignal, ph/ks(ki)-(tryph+transpose(phi))/2); legend('phase from angle - average');
% norm(ph(500:2000)/ks(ki)-(tryph(500:2000)+transpose(phi(500:2000)))/2)/norm(ph(500:2000)/ks(ki)) % for ks = 2^9
% figure; plot(collsignal(1:end-1), (ph(2:end)-ph(1:end-1))/(collsignal(2)-collsignal(1))/ks(ki)); legend('der phase from angle');
% figure; plot(collsignal, [real(signal.*exp(1i*par.k*tryph)) imag(signal.*exp(1i*par.k*tryph))]); legend('real LF
% figure; plot(collsignal(1:l), [ph(1:l)/ks(ki) tryph(1:l) transpose(phi(1:l))]);
% legend({'$\tilde{\phi}$', '$\zeta$', '$\xi$'}, 'interpreter', 'latex', 'FontSize', 15);
% xlabel('\tau_1')

figure; plot(collsignal(1:l), ph(1:l)/ks(ki), 'g'); hold on
plot(collsignal(1:l), tryph(1:l), 'b');
plot(collsignal(1:l), transpose(phi(1:l)), 'r');
legend({'$\tilde{\phi}$', '$\zeta$', '$\xi$'}, 'interpreter', 'latex', 'FontSize', 15);
xlabel('\tau_1')

%% Relative error wrt phtilde = ph/k
% symbTay = 1 + sqrt(2)*pi^2*(collsignal(1:l)-1/4).^2 -11/12*sqrt(2)*pi^4*(collsignal(1:l)-1/4).^4 ...
%     + 2783/2520/sqrt(2)*pi^6*(collsignal(1:l)-1/4).^6 -358021/205632*sqrt(2)*pi^8*(collsignal(1:l)-1/4).^8;
t1m = transpose(collsignal(1:l)-1/4);
symbTay = 1 + sqrt(2)*pi^2*t1m.^2 -11/12*sqrt(2)*pi^4*t1m.^4 ...
    + 2783/2520/sqrt(2)*pi^6*t1m.^6 -358021/205632*sqrt(2)*pi^8*t1m.^8;
% Symbolic $c_i$ & $1$  & $\sqrt{2}\pi^2$ & $\frac{-11}{12}\sqrt{2}\pi^4$ & $\frac{2783\pi^6}{2520\sqrt{2}}$ & $\frac{-358021}{205632}\sqrt{2}\pi^8$ \\

figure; plot(collsignal(1:l)-1/4, signal(1:l), 'g'); hold on
plot(collsignal(1:l)-1/4, (tryph(1:l)-ph(1:l)/ks(ki))./(ph(1:l)/ks(ki)), 'b');
plot(collsignal(1:l)-1/4, (transpose(phi(1:l))-ph(1:l)/ks(ki))./(ph(1:l)/ks(ki)), 'r');
plot(collsignal(1:l)-1/4, (symbTay-ph(1:l)/ks(ki))./(ph(1:l)/ks(ki)), 'k');
legend({'Re($V_{j,1}$)', '$(\zeta-\tilde{\phi})/\tilde{\phi}$', '$(\xi-\tilde{\phi})/\tilde{\phi}$',...
    '$(-\tilde{\phi}+\sum_{i=0}^8 c_i (\tau_1-1/4)^i)/\tilde{\phi}$'}, 'interpreter', 'latex', 'FontSize', 15);
xlabel('\tau_1')

%% Test formula for tilde{phi}
phtest = angle(signal)/ks(ki);
for i = 2:length(signal)
%     if abs(angle(signal(i) - signal(i-1) )) > 5
    if abs(angle(signal(i)) - angle(signal(i-1) ) ) > 5
%         phtest(i) = phtest(i) -2*pi*round( (phtest(i)-phtest(i-1))/2/pi);
        display(['branch cut passed at i = ' num2str(i)]);
%         phtest(i) = phtest(i) -2*pi/ks(ki)*round( angle(signal(i)-signal(i-1))/2/pi);
        phtest(i) = phtest(i) -2*pi/ks(ki)*round( angle(signal(i)) -angle(signal(i-1))/2/pi);
    end
end
phtest = phtest + 1 - phtest(closest);
phtest(closest+[-10:10])-ph(closest+[-10:10])/ks
norm(phtest-ph/ks(ki))
figure; plot(collsignal(1:l), [ph(1:l)/ks phtest(1:l)]);

%% test formula for tilde phi: Start from sp
phtest = nan*ph;
% phtest(closest) = ks(ki);
phtest(closest) = 1;
for i = 1:l/2
%     phtest(closest+i) = phtest(closest+i-1) + angle(signal(closest+i))/2/pi - angle(signal(closest+i-1))/2/pi;
%     phtest(closest+i) = phtest(closest+i-1) + angle(signal(closest+i))/2/pi/ks(ki) - angle(signal(closest+i-1))/2/pi/ks(ki);
%     phtest(closest+i) = phtest(closest+i-1) + angle(signal(closest+i)-signal(closest+i-1))/ks(ki);
%     phtest(closest+i) = phtest(closest+i-1) + (angle(signal(closest+i))-angle(signal(closest+i-1)) )/ks(ki);
    phtest(closest+i) = phtest(closest+i-1) + (angle(signal(closest+i))-angle(signal(closest+i-1)) )/ks(ki) ...
        -2*pi/ks(ki)*round( (angle(signal(closest+i)) -angle(signal(closest+i-1)))/2/pi);
%     if abs(angle(signal(closest+i)) - angle(signal(closest+i-1)) ) > 5
% %     if abs(signal(closest+i) - signal(closest+i-1)) > 5
%         display(['left branch cut passed at i = ' num2str(i)]);
% %         phtest(closest+i) = phtest(closest+i) -round( (signal(closest+i)-signal(closest+i-1))/2/pi);
%         phtest(closest+i) = phtest(closest+i) -2*pi/ks(ki)*round( (angle(signal(closest+i)) -angle(signal(closest+i-1)))/2/pi);
%     end
%     phtest(closest+i) = phtest(closest+i) -2*pi/ks(ki)*round( (angle(signal(closest+i)) -angle(signal(closest+i-1)))/2/pi);
    phtest(closest-i) = phtest(closest-i+1) + angle(signal(closest-i))/ks(ki) - angle(signal(closest-i+1))/ks(ki)...
        -2*pi/ks(ki)*floor( (angle(signal(closest-i)) -angle(signal(closest-i+1)))/2/pi);
%         -2*pi/ks(ki)*round( (angle(signal(closest-i)) -angle(signal(closest-i+1)))/2/pi);
%     if abs(angle(signal(closest-i)) - angle(signal(closest-i+1))) > 5
%         display(['right branch cut passed at i = ' num2str(i)]);
%         phtest(closest-i) = phtest(closest-i) -2*pi/ks(ki)*round( (angle(signal(closest-i)) -angle(signal(closest-i+1)))/2/pi);
%     end
end
% phtest = phtest/ks(ki);
figure; plot(collsignal(1:l), [ph(1:l)/ks phtest(1:l)]);
phtest(closest+[-10:10])-ph(closest+[-10:10])/ks
norm(phtest(1:l)-ph(1:l)/ks(ki))



%% Basic fitting
% figure; plot(collsignal(1:l) -0.25, ph(1:l)/ks(ki)); legend('phase from angle'); %-> D268 Basic fitting
% datc2 = [9.725, 12.85, 13.72, 13.911, 13.951]; figure; plot(datc2); hold on;
% datc2 = [0, -58.312, -99.94, -117.2, -122.69]; figure; plot(datc2); hold on; % Actually c4
datc2 = [-eps, -58.312, -99.94, -117.2, -122.69]; figure; plot(datc2); hold on; % Actually c4
% datc2 = [-58.312, -99.94, -117.2, -122.69]; figure; plot(datc2); hold on; % Actually c4
datc2 = [488.41, 966.1, 1231.7]; figure; plot(datc2); hold on; % Actually c6
% mat = [ones(length(datc2),1) 1./transpose(1:length(datc2))];
sgn = sign(datc2(1));
if norm(sgn-sign(datc2)) ~= 0, error('Signs not equal'); end
if 0
    deg = -3;
    mat = [ones(length(datc2),1) transpose(1:length(datc2)).^deg];
    lsq = mat\transpose(datc2)
    hold on; 
    plot(lsq(1) + lsq(2).*(1:length(datc2)).^deg)
% figure; plot(collsignal(1:l) -0.25, tryph(1:l)); legend('tau2 = 3/4');
% figure; plot(collsignal(1:l) -0.25, transpose(phi(1:l))); legend('ray perpendicular to \Gamma_2');
end
% sol = lsqnonlin(@(abc) abc(3) + exp(abc(1)-(1:length(datc2))*abc(2)) -datc2, [1,1,14])
sol = lsqnonlin(@(abc) abc(3) -sgn*exp(abc(1)-(1:length(datc2))*abc(2)) -datc2, [1,1,14])
% sol = lsqnonlin(@(abc) abc(3) - abc(1)*abc(2).^(-(1:length(datc2))) -datc2, [1,1,14])
abc = sol;
% plot(abc(3) + exp(abc(1)-(1:length(datc2))*abc(2)))
plot(abc(3) -sgn* exp(abc(1)-(1:length(datc2))*abc(2)))
% plot(abc(3) -sgn*abc(1)*abc(2).^(-(1:length(datc2) ) ) )
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
% d = 0;
d = 1;
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

%% Iterative solve
figure; subplot(1,2,1); plot(x0(1:nr)); hold on; subplot(1,2,2); plot(x0(nr+1:end)); hold on;
x = F(x0); subplot(1,2,1); plot(x(1:nr)); subplot(1,2,2); plot(x(nr+1:end));
x = F(x); subplot(1,2,1); plot(x(1:nr)); subplot(1,2,2); plot(x(nr+1:end));
x = F(x); subplot(1,2,1); plot(x(1:nr)); subplot(1,2,2); plot(x(nr+1:end));
title('\tau_2'); legend('x0','x1', 'x2', 'x3'); 
subplot(1,2,1); title('\phi(\tau_1)'); legend('x0','x1', 'x2', 'x3');


%% Load Vb1k9 etc, compute phitilde
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

ks = 2^9;
obsts = 5;

kl = length(ks);
nbOb = length(obsts);

oi = 1;
obstacle = obsts(oi);
ki = 1;
idx = (oi-1)*kl+ki;
par = getObst(obstacle); % Reset par
rx = 0.5;
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
    par.obsts(obst).hs = (par.obsts(obst).colltau(2) -par.obsts(obst).colltau(1) )/2;
end
c1x = -0.5; c1y = -0.5; c2x = -0.5; c2y = 1.5; r1 = rx; r2 = 0.5; d = 1;
nr = floor(length(par.obsts(1).colltau)/2);
collsignal = par.obsts(1).colltau(1:nr);
load Vb1k9; 
signal = Vb1k9(1:nr);

l = find(abs(collsignal-0.25) == min(abs(collsignal-0.25)));
phitilde = zeros(size(signal));
phitilde(l) = d;
for i = (l-1):-1:1
    phitilde(i) = phitilde(i+1) +(angle(signal(i))-angle(signal(i+1)))/ks(ki) -2*pi*round( (angle(signal(i))-angle(signal(i+1)))/2/pi)/ks(ki);
end
for i = (l+1):length(signal)
    phitilde(i) = phitilde(i-1) +(angle(signal(i))-angle(signal(i-1)))/ks(ki) -2*pi*round( (angle(signal(i))-angle(signal(i-1)))/2/pi)/ks(ki);
end
% figure; plot(collsignal, phitilde); xlabel('\tau_1'); ylabel('$\tilde{\phi}$', 'interpreter','latex');

%% Plot
zeta = transpose(sqrt( sum( (par.obsts(1).par(collsignal) - repmat(par.obsts(2).par(3/4), 1, length(collsignal) ) ).^2, 1) ) );
tau2Xi = (c1y-c2y+r2*sin(2*pi*collsignal) )./sqrt((c1y-c2y)^2+2*r1*sin(2*pi*collsignal)*(c1y-c2y) + r1^2);
tau2Xi = (-1).^(abs(collsignal-0.5) < 0.25).*asin(tau2Xi)/2/pi+1 -0.5*(abs(collsignal-0.5) < 0.25);
xi = sqrt( (c1x-c2x+r1*cos(2*pi*collsignal)-r1*cos(2*pi*tau2Xi) ).^2 + (c1y-c2y+r2*sin(2*pi*collsignal)-r2*sin(2*pi*tau2Xi) ).^2 );

t1m = transpose(collsignal-1/4);
symbTay0 = t1m.^0;
symbTay2 = 1 + sqrt(2)*pi^2*t1m.^2;
symbTay4 = 1 + sqrt(2)*pi^2*t1m.^2 -11/12*sqrt(2)*pi^4*t1m.^4;
% symbTay6 = 1 + sqrt(2)*pi^2*t1m.^2 -11/12*sqrt(2)*pi^4*t1m.^4 + 2783/2520/sqrt(2)*pi^6*t1m.^6;
symbTay6 = 1 + sqrt(2)*pi^2*t1m.^2 -11/12*sqrt(2)*pi^4*t1m.^4 + 2783/2520*sqrt(2)*pi^6*t1m.^6;
symbTay = 1 + sqrt(2)*pi^2*t1m.^2 -11/12*sqrt(2)*pi^4*t1m.^4 ...
    + 2783/2520*sqrt(2)*pi^6*t1m.^6 -358021/205632*sqrt(2)*pi^8*t1m.^8;
%     + 2783/2520/sqrt(2)*pi^6*t1m.^6 -358021/205632*sqrt(2)*pi^8*t1m.^8;
stot = symbTay0 + sqrt(2)*pi^2*t1m.^2.*(sqrt(2)*pi^2*t1m.^2 < 1) + -11/12*sqrt(2)*pi^4*t1m.^4.*(11/12*sqrt(2)*pi^4*t1m.^4 < sqrt(2)*pi^2*t1m.^2)...
    +2783/2520*sqrt(2)*pi^6*t1m.^6.*(2783/2520*sqrt(2)*pi^6*t1m.^6 < 11/12*sqrt(2)*pi^4*t1m.^4) ...%optimal truncation
    -358021/205632*sqrt(2)*pi^8*t1m.^8.*(358021/205632*sqrt(2)*pi^8*t1m.^8 < 2783/2520*sqrt(2)*pi^6*t1m.^6);
stot = symbTay0 + (sqrt(2)*pi^2*t1m.^2 < 1).*(sqrt(2)*pi^2*t1m.^2 + (11/12*sqrt(2)*pi^4*t1m.^4 < sqrt(2)*pi^2*t1m.^2).*(-11/12*sqrt(2)*pi^4*t1m.^4 ...
    +(2783/2520*sqrt(2)*pi^6*t1m.^6 < 11/12*sqrt(2)*pi^4*t1m.^4).*(2783/2520*sqrt(2)*pi^6*t1m.^6 ...%optimal truncation
    +(358021/205632*sqrt(2)*pi^8*t1m.^8 < 2783/2520*sqrt(2)*pi^6*t1m.^6).*(-358021/205632*sqrt(2)*pi^8*t1m.^8) ) ) );
% Symbolic $c_i$ & $1$  & $\sqrt{2}\pi^2$ & $\frac{-11}{12}\sqrt{2}\pi^4$ & $\frac{2783\pi^6}{2520\sqrt{2}}$ & $\frac{-358021}{205632}\sqrt{2}\pi^8$ \\

figure; semilogy(collsignal-1/4, (zeta-phitilde)./phitilde, 'b'); hold on;
semilogy(collsignal-1/4, (phitilde-transpose(xi))./phitilde, 'r');
semilogy(collsignal-1/4, (phitilde-symbTay0)./phitilde, 'g');
% semilogy(collsignal-1/4, (symbTay2-phitilde)./phitilde, 'm'); % Small region around SP 0.25 where negative...
semilogy(collsignal-1/4, abs(symbTay2-phitilde)./phitilde, 'm'); % Small region around SP 0.25 where negative...
semilogy(collsignal-1/4, (phitilde-symbTay4)./phitilde, 'c');
% semilogy(collsignal-1/4, (symbTay6-phitilde)./phitilde, 'r');
semilogy(collsignal-1/4, abs(symbTay6-phitilde)./phitilde, 'y');
% semilogy(collsignal-1/4, (phitilde-symbTay6)./phitilde, 'y');
semilogy(collsignal-1/4, (phitilde-symbTay)./phitilde, 'k');
% semilogy(collsignal-1/4, (phitilde-stot)./phitilde, 'r:', 'LineWidth', 2);
legend({'$(\zeta-\tilde{\phi})/\tilde{\phi}$', '$(\tilde{\phi}-\xi)/\tilde{\phi}$','$(\tilde{\phi}-c_0)/\tilde{\phi}$'...
    '$|c_0+c_2 \tau_1^2 - \tilde{\phi}|/\tilde{\phi}$', '$(\tilde{\phi}-\sum_{i=0}^4 c_i \tau_1^i)/\tilde{\phi}$', ...
     '$|\tilde{\phi}-\sum_{i=0}^4 c_i \tau_1^i|/\tilde{\phi}$', ...
    '$(\tilde{\phi}-\sum_{i=0}^8 c_i \tau_1^i)/\tilde{\phi}$'},...%,'Optimal truncation'},...%'$(\tilde{\phi}-\sum_{i=0}^T c_i \tau_1^i)/\tilde{\phi}$'}, ...
    'interpreter', 'latex', 'FontSize', 15);
%'$(c_0+c_2 \tau_1^2 - \tilde{\phi})/\tilde{\phi}$',%'$(\tilde{\phi}-\sum_{i=0}^6 c_i \tau_1^i)/\tilde{\phi}$',
xlabel('\tau_1'); ylabel('Relative error');


%% Use phi to get tau2(tau1) out of A and check with B
cl1 = transpose(collsignal);
h = cl1(2)-cl1(1);

% clos = @(phi,tau2, shft) interp1(cl1-shft*h, phi, tau2);
clos = @(phi,tau2, shft) interp1(cl1 + shft*h, phi, tau2);
dist = @(tau2) sqrt( (c1x-c2x+r1*cos(2*pi*cl1)-r1*cos(2*pi*tau2) ).^2 + (c1y-c2y+r2*sin(2*pi*cl1)-r2*sin(2*pi*tau2) ).^2 );

% F = @(x) (signal(1:nr)-dist(x)-clos(signal(1:nr), x-1/2, 0)-d );
% F = @(x) (phitilde-dist(x)-clos(phitilde, x-1/2, 0)-d );
F = @(x) (phitilde-dist(x)-clos(phitilde, x-1/2, 0)+d );
% F = @(x) ((clos(x(1:nr),x(nr+1:end)-1/2, 1)-clos(x(1:nr),x(nr+1:end)-1/2, 0))/h + ...
%     2*pi*(r1*sin(2*pi*x(nr+1:end) ).*(c1x-c2x+r1*cos(2*pi*cl1)-r1*cos(2*pi*x(nr+1:end))) ...
%     -r2*cos(2*pi*x(nr+1:end)).*(c1y-c2y+r2*sin(2*pi*cl1)+r2*sin(2*pi*x(nr+1:end))) )./dist(x(nr+1:end) ) );

x0 = 0.75*ones(size(cl1));
% x0 = 0.6+0.2*rand(1,nr);
opto = optimoptions('fsolve', 'Display', 'iter-detailed', 'FunctionTolerance', 1e-2, 'OptimalityTolerance', 1e-3, 'StepTolerance', 1e-12);
% opto = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter-detailed', 'FunctionTolerance', 1e-2*norm(F(x0)),...
%     'OptimalityTolerance', 1e-3, 'StepTolerance', 1e-3);
[XA, FVAL] = fsolve(F, x0, opto);

G = @(x) (clos(phitilde,x-1/2, 1)-clos(phitilde,x-1/2, 0))/h + (dist(x-1/2+h)-dist(x-1/2))/h ;
[norm(F(XA))/norm(F(x0)), norm(G(XA))/norm(G(x0))]

% B check A
% [XB, FVAL] = fsolve(G, x0, opto);
[XB, FVAL] = fsolve(G, XA, opto);
[norm(F(XB))/norm(F(x0)), norm(G(XB))/norm(G(x0))]
%Symbolic $a_{i+1}$ & $\frac{-1}{2\sqrt{2}+3}$ & $7\pi^2(17\sqrt{2}-24)$ & & & \\ 
%a5 == 1/84*pi^4*(1205811*sqrt(2) - 1705312), a7 == 1/128520*pi^6*(289615597399*sqrt(2) - 409578202752)
symba1 = -1/(2*sqrt(2)+3)*(cl1-1/4);
symba3 = -1/(2*sqrt(2)+3)*(cl1-1/4) +7*pi^2*(17*sqrt(2)-24)*(cl1-1/4).^3;
symba5 = -1/(2*sqrt(2)+3)*(cl1-1/4) +7*pi^2*(17*sqrt(2)-24)*(cl1-1/4).^3 +1/84*pi^4*(1205811*sqrt(2) - 1705312)*(cl1-1/4).^5;
symba7 = -1/(2*sqrt(2)+3)*(cl1-1/4) +7*pi^2*(17*sqrt(2)-24)*(cl1-1/4).^3 ...
    +pi^4/84*(1205811*sqrt(2) - 1705312)*(cl1-1/4).^5 + pi^6/128520*(289615597399*sqrt(2) - 409578202752)*(cl1-1/4).^7;
% figure; plot(cl1, [XA, XB, transpose(tau2Xi)]); ylabel('\tau_2'); xlabel('\tau_1'); legend('XA','XB', '\tau_{2,\Xi}');
figure; plot(cl1, [XA-3/4, XB-3/4, transpose(tau2Xi)-3/4, symba1, symba3, symba5, symba7]); ylabel('\tau_2'); xlabel('\tau_1'); 
legend({'$X_A$','$X_B$', '$\tau_{2,\Xi}$', '$a_1\tau_1$', '$a_1\tau_1 + a_3\tau_1^3$', '$\sum_{i=1}^5 a_i\tau_i$', ...
    '$\sum_{i=1}^7 a_i\tau_i$'}, 'interpreter', 'latex');


Xav = (XA+XB)/2;
[norm(F(Xav))/norm(F(x0)), norm(G(Xav))/norm(G(x0))]

figure; semilogy(cl1, abs([F(x0), F(XA), F(XB), G(x0), G(XA), G(XB)])); ylabel('Residu'); xlabel('\tau_1'); 
legend('F(x0)', 'F(XA)', 'F(XB)', 'G(x0)', 'G(XA)', 'G(XB)');


%% Fill series into continuous equation
Fs = [(symbTay0 -dist(3/4) - 1 +1), (symbTay2 - dist(3/4+symba1) -sqrt(2)*pi^2*symba1.^2), ...
    (symbTay4 - dist(3/4+symba3) -sqrt(2)*pi^2*symba3.^2+11/12*sqrt(2)*pi^4*symba3.^4),...
    (symbTay6 - dist(3/4+symba5) -sqrt(2)*pi^2*symba5.^2+11/12*sqrt(2)*pi^4*symba5.^4 - 2783/2520*sqrt(2)*pi^6*symba5.^6),...
    (symbTay - dist(3/4+symba7) -sqrt(2)*pi^2*symba7.^2+11/12*sqrt(2)*pi^4*symba7.^4 ...
       - 2783/2520*sqrt(2)*pi^6*symba7.^6+ 358021/205632*sqrt(2)*pi^8*symba7.^8)];
%        - 2783/2520/sqrt(2)*pi^6*symba7.^6+ 358021/205632*sqrt(2)*pi^8*symba7.^8)];
% tau_2(tau_1) -> tau2(tau-1/4)-1/2
% distc = @(tau2) sqrt( (c1x-c2x+r1*cos(2*pi*(cl1+h))-r1*cos(2*pi*tau2) ).^2 + (c1y-c2y+r2*sin(2*pi*(cl1+h))-r2*sin(2*pi*tau2) ).^2 );
% Gs = [(distc(3/4)+1-dist(3/4)-1)/h, (distc(3/4+symba1) +sqrt(2)*pi^2*(symba1+h)^2 -dist(3/4
h = 1e-7;
Gs = [(dist(3/4+h)+1-dist(3/4)-1)/h, (dist(3/4+symba1+h) +sqrt(2)*pi^2*(symba1+h).^2 -dist(3/4+symba1) -sqrt(2)*pi^2*symba1.^2)/h,....
    (dist(3/4+symba3+h) -dist(3/4+symba3) +sqrt(2)*pi^2*((symba3+h).^2-symba3.^2) -11/12*sqrt(2)*pi^4*((symba3+h).^4-symba3.^4))/h,...
    (dist(3/4+symba5+h) -dist(3/4+symba5) +sqrt(2)*pi^2*((symba5+h).^2-symba5.^2) -11/12*sqrt(2)*pi^4*((symba5+h).^4-symba5.^4) ...
       + 2783/2520*sqrt(2)*pi^6*((symba5+h).^6-symba5.^6))/h,...
    (dist(3/4+symba7+h) -dist(3/4+symba7) +sqrt(2)*pi^2*((symba7+h).^2-symba7.^2) -11/12*sqrt(2)*pi^4*((symba7+h).^4-symba7.^4) ...
       + 2783/2520*sqrt(2)*pi^6*((symba7+h).^6-symba7.^6) -358021/205632*sqrt(2)*pi^8*((symba7+h).^8-symba7.^8) )/h];
%        + 2783/2520/sqrt(2)*pi^6*((symba7+h).^6-symba7.^6) -358021/205632*sqrt(2)*pi^8*((symba7+h).^8-symba7.^8) )/h];
% symbTay = 1 + sqrt(2)*pi^2*t1m.^2 -11/12*sqrt(2)*pi^4*t1m.^4 ...
%     + 2783/2520*sqrt(2)*pi^6*t1m.^6 -358021/205632*sqrt(2)*pi^8*t1m.^8;
% figure; plot(cl1-1/4, Fs); legend('c_0', 'c_2','c_4', 'c_6','c_8'); xlabel('\tau_1'); ylabel('Error on F')
% figure; plot(cl1-1/4, Gs); legend('c_0', 'c_2','c_4', 'c_6','c_8'); xlabel('\tau_1'); ylabel('Error on G')
figure; loglog(abs(cl1-1/4), abs(Fs)); legend('c_0', 'c_2','c_4', 'c_6','c_8'); xlabel('\tau_1'); ylabel('Absolute error on F')
figure; loglog(abs(cl1-1/4), abs(Gs)); legend('c_0', 'c_2','c_4', 'c_6','c_8'); xlabel('\tau_1'); ylabel('Absolute error on G')

