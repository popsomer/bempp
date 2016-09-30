% Compute the phase of the solutions through a spectrogram and exploit this

%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

% timeRes = 32; % Higher then better localization of phase in time domain but worse in freq
timeRes = 16;
% timeRes = 64;
if 1
    ks = 2.^(7:8);
%     ks = 2^7;
%     ks = 2^8; ppw(3) = 32;
    printtoc = 7;
%     obsts = 2;
    obsts = 1;
    mti = 0;
else
    ks = 2.^(7:10);
    printtoc = 300; % Print progress info every 5min
    obsts = [3 4 8];
    mti = 2;
end
kl = length(ks);
maxob = length(obsts);
% ppw = [6 6 10 7 7 9 12 9];
ppw = zeros(8,1);
for oi = 1:8
    par = getObst(oi);
    if isfield(par, 'obsts')
        for obst = 1:length(par.obsts)
            ppw(oi) = ppw(oi) + par.obsts(obst).ppw;
        end
    else
        ppw(oi) = par.ppw;
    end
end

avm = 100; % Number of random taus to average BC over
taus = rand(avm,1);
v = struct('mti', mti, 'avm', avm, 'taus', taus, 'errBCavm', zeros(maxob*kl,2+mti),...
	'errInt', zeros(maxob*kl,2+mti), 'timeSol', zeros(maxob*kl,2+mti), 'iterGm', zeros(maxob*kl,mti), 'timeA', zeros(maxob*kl,4), 'ks', ks);

expectedEnd = 0;
for oi = 1:length(obsts)
    obstacle = obsts(oi);
	par = getObst(obstacle);
%     fN = ppw(obstacle)*256;
%     fN = ppw(obstacle)*128;
    fN = ppw(obstacle)*64;
%     fN = ppw(obstacle)*32;
    ft = linspace(0,1,fN+1); % The knots of the periodic spline;
    fco = ft(1:fN);

	for ki = 1:kl
        start = now;
		idx = (oi-1)*kl+ki;
		par.k = ks(ki);
        if isfield(par,'obsts')
            par.N = 0;
            par.r = zeros(2,length(par.obsts)); % ranges
            for obst = 1:length(par.obsts)
                par.obsts(obst).k = par.k;
                par.obsts(obst).N = ppw(obstacle)*par.k;
                par.r(1,obst) = par.N+1;
                par.N = par.N + par.obsts(obst).N;
                par.r(2,obst) = par.N;
                
                par.obsts(obst).t = linspace(0,1,par.obsts(obst).N+1); % The knots of the periodic spline;
                par.obsts(obst).colltau = par.obsts(obst).t(1:par.obsts(obst).N);
            end
        else
            par.N = ppw(obstacle)*par.k;
            par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
            par.colltau = par.t(1:par.N);
        end
        
		%% Computations
		A1 = zeros(par.N); tic; prevToc = toc;
		for i = 1:par.N
			if (toc-prevToc > printtoc)
				prevToc = toc;
				display(['Obst ' num2str(obstacle) ' at k=' num2str(ks(ki)) ': ' num2str( (i-1)/par.N,'%7.3f'),...
                    ' = part of matrix, estim sec left = ' num2str(toc*(par.N-i+1)/(i-1)) ', expected end ' datestr(expectedEnd)])
            end
            if isfield(par,'obsts')
                obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
                collx = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) );
                for obst = 1:length(par.obsts)
                    if obst == obsin
                        A1(i,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1, par.obsts(obst), 1, par.obsts(obst).N, @(taut) ones(size(taut)) );
                    else
                        A1(i,par.r(1,obst):par.r(2,obst)) = collRowQBF('error', par.obsts(obst), 1, par.obsts(obst).N, @(taut) ones(size(taut)), collx);
                    end
                end
            else
                A1(i,:) = collRowQBF(i,par,1,par.N, @(taut) ones(size(taut))); % Force not using par.difase at higher frequencies
            end
		end
		v.timeA(idx,1) = toc;
        if isfield(par,'obsts')
            b = zeros(par.N,1);
            for obst = 1:length(par.obsts)
                b(par.r(1,obst):par.r(2,obst)) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).colltau));
            end
        else
            b = par.bc(par.k,par.par(par.colltau));
        end
        c1 = A1\b; % Solution is needed for computing its phase and temporary errors
        
		%% Computing phase
		if ki == 1 % Re-use phase at higher frequencies
            tic;
            let = round(length(c1)/timeRes);
            if isfield(par, 'obsts')
                dt = 1/(par.obsts(1).colltau(2)-par.obsts(1).colltau(1));
                figure; plot([real(c1) imag(c1)])
            else
                dt = 1/(par.colltau(2)-par.colltau(1));
                figure; plot(par.colltau,[real(c1) imag(c1)])
            end
            title('Solution vector')
            [Sc,F,T] = spectrogram(c1, let,[],[], dt, 'centered');
            figure; spectrogram(c1, let,[],[], dt, 'centered'); title('Spectrogram')
%             return
            colLow = par.colltau;
            Sct = zeros(size(Sc));
            phase = zeros(size(T));
            for ti = 1:length(T)
                [~, loc] = max(Sc(:,ti));
                phase(ti) = F(loc)/par.k;
            end
            par.Tase = T;
            par.difase = arrayfun(@(Ta) quadgk(@(x) interp1(T,phase,x, 'spline'), 0, Ta), T);
			v.timeA(idx,4) = toc;
		end
		%% Compute c2
        if 0
            % Compute A2 as a standard compressed matrix
            c2 = zeros(par.N,1);
            A2 = zeros(par.N);
            tic;
            prevToc = toc;
            for i = 1:par.N
                if (toc-prevToc > printtoc)
                    prevToc = toc;
                    display(['Obst ' num2str(obstacle) ' at k=' num2str(ks(ki)) ': ' num2str(i/par.N,'%7.3f')...
                        'c2%, est. sec left=' num2str(toc*(par.N-i)/i) ', tmp solErr=' num2str(norm(c2(1:(i-1))-c1(1:(i-1)))/norm(c1(1:(i-1))  ))...
                        ', expected end ' datestr(expectedEnd)]);
                end
                tc = par.colltau(i);
                % interp1 methods: 'nearest', 'next', 'previous', 'linear','spline','pchip', or 'cubic'
                %             c2(i) = exp(1i*par.k*interp1(T,phase,tc,'spline'));
                c2(i) = exp(2i*pi*par.k*quadgk(@(x) interp1(T,phase,x,'spline'), 0, tc)); % To see whether same oscillations as c1
                A2(i,:) = collRowQBF(i,par,1,par.N);
            end % loop over row-indices i
        else
            % Compute A2 as a smaller matrix exploiting phase information
            A2 = zeros(fN);
            tN = par.N;
            par.N = fN;
            tco = par.colltau;
            par.colltau = fco;
            tt = par.t;
            par.t = ft;
            par.fco = fco;
            tic;
            prevToc = toc;
            for i = 1:fN
                tc = par.colltau(i);
                A2(i,:) = collRowQBF(i,par,1,par.N);
            end % loop over row-indices i
            par.N = tN;
            par.colltau = tco;
            par.t = tt;
        end
        
		v.timeA(idx,2) = toc;
		
		v = validate(A1,A2,par,v,idx);
% 		save('spectroPh.mat','-regexp','^(?!(A1|A2)$).')
		v.timeA(idx,3) = (now-start)*24*3600; 
        extraPol = v.timeA(idx,3)-v.timeA(idx,4); 
        powTime = 2;
        timeAllk = extraPol/(ks(ki))^powTime*sum(ks.^powTime); 
        expectedEnd = now + (timeAllk*(length(obsts)-oi) + timeAllk-extraPol/(ks(ki))^powTime*sum(ks(1:ki).^powTime) )/24/3600; 
        display(['Obstacle ' num2str(obstacle) ' at k = ' num2str(ks(ki)) ' took ' datestr(now-start, 'DD:HH:MM:SS.FFF') ...
            ' D:H:M:S and now is ' datestr(now) ', expected end ' datestr(expectedEnd)]);
	end
end
% figure; plot(par.colltau,[real(c1) 100*real(c2)])
% figure; plot(T,phase*par.k)
% c3 = A2\b; figure; plot(par.colltau, [real(c3) imag(c3)])
v
return

%% Tests FT vs phase extraction
close all
nx = 1e3;
xs = linspace(0,1,nx)';
% f = 150*(xs-0.5);
f = 1-(xs-0.5).^2;
% f = 10-5*(xs-0.5).^2;
g = (xs-0.5).^2;
k = 1e3;
v = f.*exp(1i*k*g);

figure;
let = round(nx/timeRes);
spectrogram(v, let,[],[], 1/(xs(2)-xs(1)), 'centered');
title('Spectrogram FT vs PE')
[S, F, T] = spectrogram(v, let,[],[], 1/(xs(2)-xs(1)), 'centered');
phase = zeros(size(T));
for ti = 1:length(T)
    [~, loc] = max(abs(S(:,ti)));
%     phase(ti) = F(loc)/k;
    phase(ti) = F(loc)/k*2*pi;
end
% Move the phase using knowledge of one point
[~, mi] = min(abs(xs-0.5));
% [~, mip] = min(abs(T-0.5));
% phase = phase + f(mi) - phase(mip);

% reconstr = exp(2i*pi*f.*xs);
% reconstr = exp(2i*pi*cumsum(f).*xs/k);
% reconstr = exp(2i*pi*cumsum(f)/k);
% recPh = quadgk(@(x) interp1(T,phase,x,'spline'), 0, xs);
recPh = zeros(size(xs));
reconstr = zeros(size(xs));
for xi = 1:length(xs)
%     recPh(xi) = quadgk(@(x) interp1(T,phase,x,'spline'), 0, xs(xi));
    recPh(xi) = quadgk(@(x) interp1(T,phase,x,'spline'), 0.5, xs(xi)) + g(mi); % quadgk = 0 for xs=0.5
%     reconstr(xi) = exp(2i*pi*k*quadgk(@(x) interp1(T,phase,x,'spline'), 0, xs(xi)));
%     reconstr(xi) = exp(2i*pi*k*recPh(xi) );
    reconstr(xi) = exp(1i*k*recPh(xi) );
end
reconstr = f.*reconstr;
figure; plot(xs, [g recPh]); legend('g', 'recPh');
figure; plot(xs, [real(v) real(reconstr)]); legend('v', 'reconstr');
norm(v-reconstr)/norm(v)
% norm(real(v-reconstr))/norm(real(v))
