% Compute the solution of a circle and two circles to check the extra phase.

%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

% ks = 2.^(5:9);

% bc = 1;
ks = 2^7; bc = 2;
% bcsh = pi/2; % shift in boundary condition

printtoc = 10;
kl = length(ks);
nbOb = 1;
avm = 100; % Number of random taus to average BC over
v = struct('avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(nbOb*kl,2), 'errSol', zeros(nbOb*kl,2),  ...
    'errInt', zeros(nbOb*kl,2), 'timeA', zeros(nbOb*kl,4), 'ks', ks); %, 'field', zeros(70) );

c1sin = cell(kl,1);
c1mul = cell(kl,1);

%% Computations
start = now;
for ki = 1:kl
    idx = ki;
    par = getObst(5); % Reset par
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
    
    if bc == 4
        par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/2-bcsh);sin(-pi/2-bcsh)]);
    elseif bc == 3
        par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/2);sin(-pi/2)]);
    elseif bc == 2
        par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/4);sin(-pi/4)]);
    elseif bc == 1
        par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(0);sin(0)]);
    end
    
    %% Computating full solution of multiple scattering obstacle
    A1 = zeros(par.N); tic;
    prevToc = toc;
    for i = 1:par.N
        if (toc-prevToc > printtoc)
            prevToc = toc;
            display([num2str(par.k), '=k, ' num2str( (i-1)/par.N,'%7.3f'), '=A1%, est. # sec. left for mult A1=' num2str(toc*(par.N-i+1)/(i-1)) ])
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
    b = zeros(par.N,1);
    for obst = 1:length(par.obsts)
        b(par.r(1,obst):par.r(2,obst)) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).colltau));
    end
    c1mul{ki} = A1\b;
    
    %% Computating full solution of single scattering obstacle
    par = getObst(1); % Reset par
    par.ppw = 10;
    par.k = ks(ki);
    par.N = par.ppw*par.k;
    par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
    par.colltau = par.t(1:par.N);
    
    if bc == 4
        par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/2-bcsh);sin(-pi/2-bcsh)]);
    elseif bc == 3
        par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/2);sin(-pi/2)]);
    elseif bc == 2
        par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/4);sin(-pi/4)]);
    elseif bc == 1
        par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(0);sin(0)]);
    end
    
    A1 = zeros(par.N); tic;
    prevToc = toc;
    for i = 1:par.N
        % 	parfor i=1:par.N % Instead of a sequential loop
        if (toc-prevToc > printtoc)
            prevToc = toc;
            display([num2str(par.k), '=k, ' num2str( (i-1)/par.N,'%7.3f'), '=A1%, est. # sec. left for single A1=' num2str(toc*(par.N-i+1)/(i-1)) ])
        end
        A1(i,:) = collRowQBF(i,par);
    end
    b = par.bc(par.k,par.par(par.colltau));
    c1sin{ki} = A1\b; % The solution is needed for computing the correlations.
    
%     save c1OneCirc.mat
    display([num2str(ki) '=ki, now is ' datestr(now) ', expected end ' datestr(start + ...
        (now-start)*sum(ks.^2)/sum(ks(1:ki).^2) )  ]);
end


return


%% Substract solution from one circle with phase shift
load c1OneCirc.mat

%% ki
ki = 4;

c1m = c1mul{ki};
c1s = c1sin{ki};
len = length(c1s);
slice = 1:len;

c1slice = c1m(slice);
% figure; plot([real(c1slice) imag(c1slice)]); title('c1slice');
% figure; plot([real(c1s) imag(c1s)]); title('c1single');

% figure; plot([real(c1slice./c1s) imag(c1slice./c1s)]); title('c1slice./c1s');
if bc == 1
    rat = mean(c1slice(round(60/128*len):round(115/128*len))./c1s(round(60/128*len):round(115/128*len)))
elseif bc == 2
    rat = mean(c1slice(round(60/128*len):round(90/128*len))./c1s(round(60/128*len):round(90/128*len)))
else
    rat = nan;
end
c1sp = c1slice - c1s*rat;
figure; plot([real(c1sp) imag(c1sp)]); title('c1sp');
% figure; spectrogram(c1slice, round(length(c1sp)/16),[],[], 1/len, 'centered'); title('c1slice')
% figure; spectrogram(c1sp, round(length(c1sp)/16),[],[], 1/len, 'centered');title('c1sp')

par = getObst(5);
N = par.obsts(1).ppw*ks(ki);
colltau = linspace(0,1,N+1);
colltau = colltau(1:N);

%% Extract phase from angle
% collsignal = colltau(1:(length(colltau)/4));
collsignal = colltau(1:(length(colltau)/2));
signal = c1sp(1:(length(colltau)/2));
figure; plot(collsignal, [real(signal) imag(signal)]);
l = length(signal)/2;
% ph = nan*signal;
ph = zeros(size(signal));
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
figure; plot(collsignal, ph/ks(ki));

%% Guess phase
hold on;
trypt = [par.obsts(2).par(3/4); 1.05];
trypt(1) = trypt(1) - 0.31;
plot(collsignal,sqrt( sum( (par.obsts(1).par(collsignal) - repmat(trypt(1:2), 1, length(collsignal) ) ).^2, 1) )-trypt(3));
% plot(collsignal,sqrt( sum( (par.obsts(1).par(collsignal) - repmat(par.obsts(2).par(3/4), 1, length(collsignal) ) ).^2, 1) )-1);
% plot(collsignal,sqrt( sum( (par.obsts(1).par(collsignal) - repmat(par.obsts(1).par(1/4)/2+par.obsts(2).par(3/4)/2, 1, length(collsignal) ) ).^2, 1) )-0.5);
% ptshsol = fminunc(@(ptsh) norm(transpose(ph/ks(ki)) - sqrt( sum( (par.obsts(1).par(collsignal) - repmat(ptsh(1:2), 1, length(collsignal) ) ).^2, 1) ) - ptsh(3)), ...
%     trypt);
trypt = [par.obsts(2).par(3/4); 1.05; 1/2; -1/2; 1; 1];
[ptshsol, fcv] = fminunc(@(ptsh) norm(transpose(ph/ks(ki)) - (ptsh(7)*sqrt( sum( (par.obsts(1).par(collsignal) - repmat(ptsh(1:2), 1, length(collsignal) ) ).^2, 1) ) ...
    - ptsh(3) + ptsh(4)*cos(2*pi*collsignal) + transpose(ptsh(5)*cos(2*pi*tau2s) + ptsh(6)*posphas)) ), trypt);
%     [par.obsts(2).par(3/4); 1] );
plot(collsignal,sqrt( sum( (par.obsts(1).par(collsignal) - repmat(ptshsol(1:2), 1, length(collsignal) ) ).^2, 1) ) + ptshsol(3));

%% Compute phase from d(x_2+d)/dtau_2 = 0
tau1 = 1/3;
t2 = fminbnd( @(tau2) cos(2*pi*tau2)/2 + sqrt(sum( (par.obsts(1).par(tau1) - par.obsts(2).par(tau2) ).^2, 1) ), 0.5, 1)
posphas = 0*ph;
% tau2s = 0*ph;
tau1s = 0*ph;
coll2 = colltau((length(colltau)/2 + 1):end);
for i=1:length(signal)
%     [t2, fval] = fminbnd( @(tau2) cos(2*pi*tau2)/2 + sqrt(sum( (par.obsts(1).par(collsignal(i) ) - par.obsts(2).par(tau2) ).^2, 1) ), 0.5, 1);
%     [t2, fval] = fminbnd( @(tau2) sqrt(sum( (par.obsts(1).par(collsignal(i) ) - par.obsts(2).par(tau2) ).^2, 1) ), 0.5, 1);
%     [t1, fval] = fminbnd( @(tau1) sqrt(sum( (par.obsts(1).par(tau1) - par.obsts(2).par(coll2(i)) ).^2, 1) ), 0, 0.5); % Should be SP/refl on obst 1
    [t1, fval] = fminbnd( @(tau1) cos(2*pi*tau1)/2 + sqrt(sum( (par.obsts(1).par(tau1) - par.obsts(2).par(coll2(i)) ).^2, 1) ), 0, 0.5);
%     posphas(i) = fval;
%     posphas(i) = cos(2*pi*coll2(i))/2- sqrt(sum( (par.obsts(1).par(t1) - par.obsts(2).par(coll2(i)) ).^2, 1) );
    posphas(i) = sqrt(sum( (par.obsts(1).par(t1) - par.obsts(2).par(coll2(i)) ).^2, 1) );
%     posphas(i) = fval-cos(2*pi*t2)/2;
%     posphas(i) = fval-cos(2*pi*collsignal(i))/2;
%     tau2s(i) = t2;
    tau1s(i) = t1;
end
figure; plot(coll2, tau1s);
% figure; plot(collsignal, tau2s);
% figure; plot(collsignal, posphas);
% figure; plot(collsignal, ph/ks(ki)); hold on; plot(collsignal, posphas-0.8);

%% Use phase
% for i = 1:length(signal)
%     posphas(i) = posphas(i) - sqrt(sum( (par.obsts(1).par(collsignal(i) ) - par.obsts(2).par(tau2s(i)) ).^2, 1) );
%     posphas(i) = posphas(i) + cos(2*pi*collsignal(i))/2;
%     posphas(i) = posphas(i) - cos(2*pi*tau2s(i))/2;
%     posphas(i) = posphas(i) - cos(2*pi*tau1s(i))/2;
%     posphas(i) = posphas(i) + cos(2*pi*coll2(i))/2;
%     posphas(i) = posphas(i) - sqrt(sum( (par.obsts(1).par(tau1s(i)) - par.obsts(2).par(coll2(i)) ).^2, 1) );
% end


if 0
    phas = exp(1i*ks(ki)*sqrt( sum( (par.obsts(1).par(colltau) - repmat(par.obsts(2).par(3/4), 1, length(colltau) ) ).^2, 1) ));
    tmp = c1sp.*transpose(1./phas);
    figure; plot([real(tmp) imag(tmp)]); title([num2str(ks(ki)) ' = k, tmp, ki = ' num2str(ki)]);
    % figure; plot(abs(fft(tmp))); title('ffttmp');
elseif 1
%     phas = transpose(exp(1i*ks(ki)*ph));
    phas = transpose(exp(1i*ph));
    tmp = c1sp(1:length(ph)).*transpose(1./phas);
    figure; plot(collsignal, [real(tmp) imag(tmp)]); title([num2str(ks(ki)) ' = k, tmp from ph, ki = ' num2str(ki)]);
else
%     phas = transpose(exp(1i*ks(ki)*(posphas-0.85) ));
    phas = transpose(exp(1i*ks(ki)*posphas));
    tmp = c1sp(1:length(posphas)).*transpose(1./phas);
    figure; plot(collsignal, [real(tmp) imag(tmp)]); title([num2str(ks(ki)) ' = k, tmp, ki = ' num2str(ki)]);
end

%% Compute extra phase
signal = tmp(1:length(colltau)/2);
collsignal = colltau(1:(length(colltau)/2));
l = length(signal)/2;

extraph = zeros(size(signal));
for i = l:-1:1
    extraph(i) = angle(signal(i));
    if (i ~= l) && (abs(extraph(i) - extraph(i+1)) > 5)
        extraph(i) = extraph(i) -2*pi*round( (extraph(i)-extraph(i+1))/2/pi);
    end
end
for i = l:length(signal)
    extraph(i) = angle(signal(i));
    if (i ~= l) && (abs(extraph(i) - extraph(i-1)) > 5)
        extraph(i) = extraph(i) -2*pi*round( (extraph(i)-extraph(i-1))/2/pi);
    end
end
figure; plot(collsignal, extraph/ks(ki));
