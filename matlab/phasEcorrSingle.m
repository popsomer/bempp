% Compute the phase from the correlations for single scattering obstacles.

if 1
    
%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

percDecay = 1; % Percentage of the window for the C-inf decay: 0 means a block window and 1 means not identically one on any interval.

thrType = 'l'; % Or 'g' for global threshold
Tcor = 0.02;
% If corrDist is nonzero, use the correlations iso the physical distance to determine the window. The value gives 
% the percentage below the threshold from where the window is identically zero.
corrDist = 0; 

% ks = 2.^7;
% ks = 2.^(7:10);
ks = 2^7;
obsts = 1;
% obsts = 9;
bc = 1;

mti = 0;
printtoc = 3;
kl = length(ks);
nbOb = length(obsts);

avm = 100; % Number of random taus to average BC over
v = struct('mti', mti, 'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(nbOb*kl,2+mti),...
	'perc', zeros(nbOb*kl,2), 'errSol', zeros(nbOb*kl,2+mti),  'errInt', zeros(nbOb*kl,2+mti), ...
	'timeSol', zeros(nbOb*kl,2+mti), 'nbIter', zeros(nbOb*kl,mti), 'timeA', zeros(nbOb*kl,4), 'ks', ks);

%% Computations
for oi = 1:length(obsts)
	obstacle = obsts(oi);
	par = getObst(obstacle);   
        
	start = now;
%     fN = par.ppw*20;
%     fN = par.ppw*64;
%     fN = par.ppw*ks(1);
%     ft = linspace(0,1,fN+1); % The knots of the periodic spline;
%     fco = ft(1:fN);
    
	for ki = 1:kl
		idx = (oi-1)*kl+ki;
        par = getObst(obstacle); % Reset par
%     par.ppw = 10;
		par.k = ks(ki);
		par.N = par.ppw*par.k;
        par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
        par.colltau = par.t(1:par.N);
        if bc == 2
            par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/4);sin(-pi/4)]);
        elseif bc == 1
            par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(0);sin(0)]);
        end
        
        if isfield(par,'phase')
            par = rmfield(par, 'phase'); % no phase when computing A1
%             par = rmfield(par, 'quadR'); % don't need additional points when computing A1
        end
        qr1 = 12;
%         par.quadR = qr1; % Seems to give oscillations of c1 in the shadow region
        
		%% Computating full solution
		A1 = zeros(par.N); tic;
        prevToc = toc;
        for i = 1:par.N
            % 	parfor i=1:par.N % Instead of a sequential loop
            if (toc-prevToc > printtoc)
                prevToc = toc;
                display([num2str(par.k), '=k, ' num2str( (i-1)/par.N,'%7.3f'), '=A1%, est. # sec. left for A1=' num2str(toc*(par.N-i+1)/(i-1)) ])
            end
			A1(i,:) = collRowQBF(i,par);
		end
        v.timeA(idx,1) = toc;
        b = par.bc(par.k,par.par(par.colltau));
        c1 = A1\b; % The solution is needed for computing the correlations.
        
        %% Computing correlations
%         if ki == 1 % Re-use the R that is computed here at higher frequencies.
        if 0 % Re-use the R that is computed here at higher frequencies.
            tic
%             [R, sigma] = calcCorr(par, c1, Tcor, percDecay); % Don't use A1, but the integral.
            [R, sigma] = calcCorr(par, c1, Tcor, percDecay, [], A1);
            v.timeA(idx,4) = toc;
            colLow = par.colltau; % Save these for higher freqencies.
        end
		%% Compute A2
%         par.phase = @(t) transpose(par.par(t) )*[cos(0);sin(0)]; % For incident plane wave in x-direction; only set phase after computing A1
        if bc == 1 %obstacle <= 1
            par.phase = @(t) transpose(transpose(par.par(t) )*[cos(0);sin(0)]); % For incident plane wave in x-direction; only set phase after computing A1
        elseif bc == 2 %obstacle == 2
            par.phase = @(t) transpose(transpose(par.par(t) )*[cos(-pi/4); sin(-pi/4)]); % For incident plane wave in x-direction; only set phase after computing A1
        end
%         A2 = nan*zeros(par.N);
        if 0
            A2 = zeros(par.N); % We could actually use a sparse matrix structure now because its structure is given by R.
            curThr = par.xi*max(max(abs(R)));
            tic;
            for i = 1:par.N
                tc = par.colltau(i);
                [~, cli] = min(abs(colLow-tc)); % Index i in R closest to tc
                if thrType == 'l', curThr = par.xi*max(abs(R(cli,:))); end
                if corrDist
                    toAdd = Cinf(sigma,tc-Tcor, tc-Tcor/5, tc+Tcor/5, tc+Tcor,1); % Add to the diagonal to ensure the Green singularity is included.
                    A2(i,:) = windRow(i, par, [], [], 0, 1, abs(R(i,:)/curThr) +toAdd, corrDist, sigma);
                else % Physical distance
                    I = find(abs(R(cli,:)) >= curThr);
                    ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
                    bounds = [tc-Tcor, sigma(I(ini(1:(length(ini)-1) )+1)); tc+Tcor, sigma(I(ini(2:length(ini))))]; % identically 1 on tc \pm T
                    decay = repmat(Tcor*percDecay,1,size(bounds,2));
                    bounds = [(bounds + [-decay; decay]); [1;1]*decay];
                    A2(i,:) = windRow(i,par,bounds);
                end
            end % loop over row-indices i
            v.timeA(idx,2) = toc;
        elseif 1
            fact = 2^(1/4);
            %             for iter = 1:log(par.k/8)/log(fact)
            %                 fN = round(8*par.ppw*fact^iter);
            for iter = 1:1%log(par.k/24)/log(fact)
%                 fN = round(24*par.ppw*fact^iter);

                %% Setting up compression
                fN = 48; %round(par.N/32);
                ft = linspace(0,1,fN+1); % The knots of the periodic spline;
                fco = ft(1:fN);
                
%                 par.quadR = max(2, round(par.N/fN) );
%                 par.quadR = max(2, par.quadR*round(par.N/fN) );
%                 par.quadR = max(2, qr1*round(par.N/fN) );
                par.quadR = max(2, 6*round(par.N/fN) );
%                 par.quadR = 20;
                
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
                    if (toc-prevToc > printtoc)
                        prevToc = toc;
                        display([num2str(par.k), ' = k, ' num2str( (i-1)/fN,'%7.3f'), ' = A2%, est. # sec. left for A2 = ' num2str(toc*(fN-i+1)/(i-1)) ])
                    end
                    tc = par.colltau(i);
                    %                 A2(i,:) = collRowQBF(i,par,1,par.N);
                    A2(i,:) = collRowQBF(i,par);
                end % loop over row-indices i
                par.N = tN;
                par.colltau = tco;
                par.t = tt;
                
                b2 = par.bc(par.k,par.par(par.fco));
                %                 c2 = A2\b2;
                if ki ~= 1
%                     c2 = ks(ki)/ks(ki-1)*interp1([oldfco (oldfco(1) + 1)], [oldc2; oldc2(1)], fco'); % Interpolate old c2
%                     c2 = gmres(A2,b2,fN, 1e-5, fN, [], [], c2);
%                     c2 = gmres(A2,b2,fN, 1e-5, fN);
                    c2 = A2\b2;
                else
                    c2 = A2\b2;
                end
%                 c2 = +b2*par.k*2i.*cos(2*pi*fco').*(fco' > 0.25).*(fco' < 0.75);
%                 c2 = -par.k*2i.*cos(2*pi*fco').*(fco' > 0.25).*(fco' < 0.75);
                erSol = norm(interp1([par.fco (par.fco(1)+1)], [c2; c2(1)], par.colltau').*transpose(exp(1i*par.k*par.phase(par.colltau) ))-c1)/norm(c1);
                
                %% finishing
%                 figure; plot(par.fco, [real(c2) imag(c2)]); legend('Re(c2) using fco','Im(c2) using fco')
                figure; plot(par.colltau, par.phase(par.colltau)); title('Phase')
                figure; spectrogram(c1, round(length(c1)/16),[],[], 1/(par.colltau(2)-par.colltau(1)), 'centered'); title('Spectrogram')
                figure; plot(par.colltau, [real(c1) imag(c1)]); legend('Re(c1) using fco','Im(c1) using fco')
                %                 v = validate(A1,A2,par,v,idx);
                %                 erSol = v.errSol(ki,2);
                display([num2str(ks(ki)) ' = k, iter = ' num2str(iter) ', fN = ' num2str(fN) ' and erSol = ' num2str(erSol) ])
                %                 display([num2str(ks(ki)) ' = k, iter = ' num2str(iter) ', fN = ' num2str(fN) ...
                %                     ', errBCavm = ' num2str(v.errBCavm(ki,:) ) ' and erSol = ' num2str(erSol) ])
                if erSol < 0.04; %0.01
                    break;
                end
            end
        end
        oldfco = fco;
        oldc2 = c2;
        if ki == 1
            firstfco = fco;
            firstc2 = c2;
        elseif ki == 2
            secondfco = fco;
            secondc2 = c2;
        end
        v = validate(A1,A2,par,v,idx)
        display([num2str(oi) ' = oi, ki = ' num2str(ki) ', now is ' datestr(now) ', expected end ' datestr(start + ...
            (now-start)*sum(ks.^2)/sum(ks(1:ki).^2)*( length(obsts)-oi + 1) )  ]);
	end
end
% v
return

if 0
    figure;
    lambdas = logspace(-5,5,10);
    legs = cell(length(lambdas)+2,1);
%     plot(fco, [real(c2fromfirst) real(A2\b2)]);
    plot(fco, imag([c2fromfirst A2\b2]));
    hold on;
    legs{1} = 'c_f';
    legs{2} = 'A2 \\ b2';
    for li = 1:length(lambdas)
        c2l = (lambdas(li)*transpose(A2)*A2 + eye(fN) )\(lambdas(li)*transpose(A2)*b2 + c2fromfirst);
%         plot(fco, real(c2l) );
        plot(fco, imag(c2l) );
%         hold on;
        legs{li+2} = ['\lambda = ' num2str(lambdas(li))];
    end
    legend(legs)
%     legend(lambdas');
% y = fft(c2fromFine); f = (0:length(y)-1)/(fco(2)-fco(1))/length(y); figure; subplot(2,1,1); plot(f,abs(y) ); title('Magnitude'); subplot(2,1,2); plot(f,angle(y)); title('Phase')

end

%% Show the correlations
if 0
figure; surf(abs(R), 'EdgeColor','none'); xlabel('n'); ylabel('m'); 

fs = 20;
figure; 
% pcolor(min(abs(R),0.1) );
pcolor(abs(R) );
hold on
shading interp;
xlabel('n'); 
ylabel('m'); 
set(gca,'FontSize',fs);

hh = colorbar(); 
% ylabel(hh,'min$(|r_{m,n}|,0.1)$', 'interpreter','latex', 'FontSize',fs); 
ylabel(hh,'$|r_{m,n}|$', 'interpreter','latex', 'FontSize',fs); 
set(gca,'YDir','rev')

figure; surf(angle(R), 'EdgeColor','none'); xlabel('n'); ylabel('m');

figure; plot(sigma, [real(R(round(size(R,1)*0.5), :)); imag(R(round(size(R,1)*0.5), :))] );
end

end


%% R around diagonal (or nonzero contrib)
idSig = 0*par.colltau;
Rpl = idSig;
dist = Rpl;
is = [round(par.N/2):par.N, 1:(round(par.N/2)-1)];
phIll = Rpl;
phase = Rpl;
% for i = 1:par.N
for nr = 1:length(is)
    i = is(nr);
%     [~, ix] = min( (par.colltau(i)-sigma).^2);
    [~, ix] = max(abs(R(i,:)));
    idSig(i) = ix;
    Rpl(i) = R(i,ix);
%     Rpl(i) = b(i);
%     Rpl(i) = c1(i);
    dist(i) = norm(par.par(par.colltau(i)) - par.par(sigma(ix)));
    phIll(i) = par.k*transpose(par.par(par.colltau(i)))*[cos(0);sin(0)];
    phase(i) = angle(Rpl(i));
    if nr == 1
        display([num2str(phase(i)) ' = phase in illum and phase of BC = ' num2str(angle(b(i)))]);
    elseif phase(i) - phase(is(nr-1)) > 5
%         phase(i) = phase(i) -2*pi;
        phase(i) = phase(i) -2*pi*round( (phase(i)-phase(is(nr-1)))/2/pi);
    elseif phase(is(nr-1)) - phase(i) > 5
%         phase(i) = phase(i) +2*pi;
%         phase(i) = phase(i) +2*pi*ceil( (phase(is(nr-1))-phase(i))/2/pi);
        phase(i) = phase(i) +2*pi*round( (phase(is(nr-1))-phase(i))/2/pi);
    end
end
% figure; plot(abs(Rpl));
% figure; plot(angle(Rpl)-dist);
figure; plot(transpose([phase; phIll]));

%% Short-time averages
rad = 38;
rad = 58;
dtau = 2*par.colltau(rad);
c3 = nan*c1;
for i = (rad+1):par.N-rad
    c3(i) = sum(exp(-1i*par.k*par.phase(par.colltau((i-rad):(i+rad))))*c1((i-rad):(i+rad)))/(2*rad+1);
end
figure; plot(par.colltau, [real(c3) imag(c3)]); title('c3');
c4 = c1-c3.*transpose(exp(-1i*par.k*par.phase(par.colltau)));
figure; plot(par.colltau, [real(c4) imag(c4)]); title('c4');

tesph = @(t) sqrt(sum( (par.par(t) - par.par(0.4-t)).^2,1));
c5 = nan*c1;
for i = (rad+1):par.N-rad
%     c5(i) = sum(exp(-1i*par.k*tesph(par.colltau((i-rad):(i+rad))))*c1((i-rad):(i+rad)))/(2*rad+1);
    c5(i) = sum(exp(-1i*par.k*tesph(par.colltau((i-rad):(i+rad))))*c4((i-rad):(i+rad)))/(2*rad+1);
end
figure; plot(par.colltau, [real(c5) imag(c5)]); title('c5');

%% Fourier transform removal HF
% length(fft)/2 = alternating +-1 and we have ppw samples per period
% cut = ceil(par.ppw/4);
cut = ceil(length(c1)/par.ppw/4);
tmp = c1.*transpose(exp(-1i*par.k*par.phase(par.colltau)));
ffttmp = fft(tmp);
ffttmp((cut+1):(end-cut)) = 0;
c3 = ifft(ffttmp);
figure; plot(par.colltau, [real(c3) imag(c3)]); title('c3');
% c4 = c1-c3.*transpose(exp(-1i*par.k*par.phase(par.colltau)));
c4 = c1-c3.*transpose(exp(1i*par.k*par.phase(par.colltau)));
figure; plot(par.colltau, [real(c4) imag(c4)]); title('c4');
figure; spectrogram(c4, round(length(c4)/16),[],[], 1/(par.colltau(2)-par.colltau(1)), 'centered');

%% Find phase by 'ray tracing'
sp = zeros(par.N,2);
ph = zeros(par.N,2);
z = par.par(linspace(0,1,100)); figure; plot(z(1,:), z(2,:))
for i = round(par.N*0.2):round(par.N*0.4)
%     sp(i) = 0.4;
    spi = round(par.N*0.4);
    phix = 1;
    sgn = 1;
    while (phix < 3) && (spi < round(par.N*0.6) )
        dif = par.par(par.colltau(spi)) - par.par(par.colltau(i));
        alpha = atan2(dif(2), dif(1) );
        norma = par.normal(par.colltau(spi));
        beta = atan2(norma(2), norma(1));
%         delta = beta-pi+alpha
        if sgn*(mod(2*beta-pi+alpha, 2*pi) - pi) > 0
            sp(i, phix) = par.colltau(spi);
%             ph(i, phix) = norm(par.colltau(i)-par.colltau(spi));
%             ph(i, phix) = norm(par.colltau(i)-par.colltau(spi)) + [1, 0]*par.par(par.colltau(spi)); % Add phase of incident wave at sp
            ph(i, phix) = norm(par.colltau(i)-par.colltau(spi)) - [1, 0]*par.par(par.colltau(spi)); % Add phase of incident wave at sp
            phix = phix + 1;
            sgn = -sgn;
        end
        spi = spi + 1;
    end
end
for i = round(par.N*0.4):round(par.N*0.6)
    spi = round(par.N*0.4);
    phix = 1;
    sgn = -1;
    while (phix < 3) && (spi > round(par.N*0.2) )
        dif = par.par(par.colltau(spi)) - par.par(par.colltau(i));
        alpha = atan2(dif(2), dif(1) );
        norma = par.normal(par.colltau(spi));
        beta = atan2(norma(2), norma(1));
%         delta = beta-pi+alpha
        if sgn*(mod(2*beta-pi+alpha, 2*pi) - pi) > 0
            sp(i, phix) = par.colltau(spi);
%             ph(i, phix) = norm(par.colltau(i)-par.colltau(spi));
%             ph(i, phix) = norm(par.colltau(i)-par.colltau(spi)) + [1, 0]*par.par(par.colltau(spi)); % Add phase of incident wave at sp
            ph(i, phix) = norm(par.colltau(i)-par.colltau(spi)) - [1, 0]*par.par(par.colltau(spi)); % Add phase of incident wave at sp
            phix = phix + 1;
            sgn = -sgn;
        end
        spi = spi - 1;
    end
end
figure; plot(par.colltau, sp);

%% remove from already removed inc phase
ixph = 2;
sgn = 1;
tmp = c4.*exp(-sgn*1i*par.k*ph(:,ixph));
ffttmp = fft(tmp);
ffttmp((cut+1):(end-cut)) = 0;
c5 = ifft(ffttmp);
figure; plot(par.colltau, [real(c5) imag(c5)]); title('c5');
c6 = c4-c5.*exp(sgn*1i*par.k*ph(:,ixph) );
figure; plot(par.colltau, [real(c6) imag(c6)]); title('c6');
figure; spectrogram(c6, round(length(c6)/16),[],[], 1/(par.colltau(2)-par.colltau(1)), 'centered');

