% Introduce window functions with a wavenumber-dependent support to reduce the complexity of a standard BEM for a circle.
%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

percDecay = 0.5; % Percentage of the window for the C-inf decay: 0 means block window and 1 means not identically one on any interval
ks = 2.^(4:12);

kl = length(ks);
mti = 2; % gmres on A1\b, A2\b, A1\b with precond A2, A1\b with x0=A2\b

avm = 100; % Number of random taus to average the error on the BC over
v = struct('conds', zeros(kl,2), 'mti', mti, 'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(kl,2+mti),...
    'perc', zeros(kl,2), 'errSol', zeros(kl,2+mti), 'compresErr', zeros(kl,2), 'timeSol', zeros(kl,2+mti), ...
    'nbIter', zeros(kl,mti), 'timeA', zeros(kl,2), 'ks', ks, 'errInt', zeros(kl,2+mti) ); 

%% Computations
for oi = 0:1 % Cubic and linear basis functions
    par = getObst(oi);
    start = now;
    for ki = 1:kl
        par.k = ks(ki); par.N = round(par.ppw*par.k);
        par.t = linspace(0,1,par.N+1); % The knots of the periodic spline;
        par.colltau = par.t(1:par.N);
        
        A1 = zeros(par.N);
        tic;
        for i=1:par.N % Instead of a sequential loop
            A1(i,:) = collRowQBF(i,par);
        end
        v.timeA(ki,1) = toc;
        tic;
        
        fact = 0.4;
        TtrRow = fact*0.7*(par.k/64)^(-1/3);
        TtrColu = fact*1.2*(par.k/64)^(-1/3);
        Tsp = fact*1.7*(par.k/64)^(-1/2);
        Tsi = fact*0.35*(par.k/64)^(-1);
        skew = 13;
        A2 = zeros(par.N);
        
        tic;
        for i = 1:par.N
            tc = par.colltau(i);
            if (abs(tc-0.25) < TtrRow)
                tbounds = [tc-Tsi; max(tc+TtrColu, 0.5-tc+TtrColu); Tsi*percDecay; TtrColu*percDecay]; % Much faster oscillations to the left
            elseif (abs(tc-0.75) < TtrRow)
                tbounds = [min(tc-TtrColu, 1.5-tc-TtrColu); tc+Tsi; TtrColu*percDecay; Tsi*percDecay];
            elseif (tc > 0.25) && (tc < 0.75)
                tbounds = [tc-Tsi-Tsi*skew*(tc-0.25)/0.5; tc+Tsi+Tsi*skew*(0.75-tc)/0.5; ...
                    percDecay*(Tsi+Tsi*skew*(tc-0.25)/0.5); percDecay*(Tsi+Tsi*skew*(0.75-tc)/0.5)]; % 10*Tsi to the right of tc when at/near 0.25
            elseif (tc < 0.25)
                tbounds = [tc-Tsi, 0.5-tc-Tsp; tc+Tsi, 0.5-tc+Tsp;...
                    percDecay*Tsi, percDecay*Tsp; percDecay*Tsi, percDecay*Tsp]; % Shadow: singularity and stationary point
            else
                tbounds = [1.5-tc-Tsp, tc-Tsi; 1.5-tc+Tsp, tc+Tsi;...
                    percDecay*Tsp, percDecay*Tsi; percDecay*Tsp, percDecay*Tsi]; % Shadow: singularity and stationary point
            end
            tbounds(1,:) = tbounds(1,:) - eps;
            tbounds(2,:) = tbounds(2,:) + eps;
            A2(i,:) = windRow(i,par,tbounds);
        end % loop over row-indices i
        v.timeA(ki,2) = toc;
        v = validate(A1,A2,par,v,ki);
        
        clearvars A1 A2 % Avoid saving large matrices
        save complRed.mat
        if oi == 0
            vCub = v;
        end
        % Assume constructing the matrix dominates, so O(k^2), and linear basis functions are about 11 times faster.
        display([num2str(oi) ' = oi, ki = ' num2str(ki) ', now is ' datestr(now) ', expected end ' datestr(start + ...
            (now-start)*sum(ks.^2)/sum(ks(1:ki).^2)*( (1-oi)/11 + 1) )  ]);
    end % ki loop
end % degree basis functions loop


%% Make the plots shown in the article for cubic and linear basis functions.
lws = 'LineWidth'; lw = 5;
fss = 'Fontsize'; fs = 22;
mss = 'MarkerSize'; ms = 15;
l = {'v', '+', 'o', 'x', '*', 'h', 'd'};
ll = {'-','--',':', '-.'};
c = 'bgrkcmy';
figure;
loglog(v.ks, v.conds(:,1), [l{1} c(1)], mss,ms); hold on;
loglog(v.ks, v.conds(:,2), [l{2} c(2)], mss,ms); hold on;

loglog(vCub.ks, vCub.conds(:,1), [l{3} c(3)], mss,ms); hold on;
loglog(vCub.ks, vCub.conds(:,2), [l{4} c(4)], mss,ms); hold on;
xlabel('k',fss,fs);
legend({'$A$ (L)', '$\tilde{A}$ (L)', '$A$ (C)', '$\tilde{A}$ (C)'}, 'interpreter','latex');
ylabel('Condition number',fss,fs);
set(gca,fss,fs);

figure;
loglog(v.ks, v.nbIter(:,1), [l{1} c(1)], mss,ms); hold on;
loglog(v.ks, v.nbIter(:,2), [l{2} c(2)], mss,ms); hold on;
loglog(vCub.ks, vCub.nbIter(:,1), [l{3} c(3)], mss,ms); hold on;
loglog(vCub.ks, vCub.nbIter(:,2), [l{4} c(4)], mss,ms); hold on;
xlabel('k',fss,fs);
legend({'$A$\textbackslash$b$ (L)', '$\tilde{A}$\textbackslash$b$ (L)', ...
    '$A$\textbackslash$b$ (C)', '$\tilde{A}$\textbackslash$b$ (C)'}, 'interpreter','latex');
ylabel('# Gmres iterations',fss,fs);
set(gca,fss,fs);
