% Introduce window functions with fixed support in a standard BEM for a circle.
%% Initialising
clear variables -except A3 A
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

% If varykTs = 0 then plot sparsity structure for one k and T, else vary k and T
varykT = 0;
factShad = 3; % Factor of the width of a window around a stationary point, although this is dependent on k.
percDecay = 0.5; % Percentage of the window for the C-inf decay: 0 means a block window and 1 means not identically one on any interval.

par = getObst(1);

if varykT
    ks = 2.^(4:10); Ts = linspace(0.001,0.15,20); mti = 2;
%     ks = 2.^(4:6); Ts = linspace(0.001,0.15,5);
else % Only one k and T
    ks = 2^8; Ts = 0.08; mti = 0;
%     ks = 2^4; Ts = 0.001; mti = 0;
end
Tl = length(Ts); kl = length(ks);
avm = 100; % Number of random taus to average BC over
taus = rand(avm,1);
rtests = 1 + 2*(randn(avm,1)).^2; % To compare with exact solution in the field

% v is for validation: 2 refers to A and A2=\tilde{A}, 2+mti to A\b, A2\b and the solutions of the iterative solver
if varykT
    v = struct('conds', zeros(Tl*kl,2), 'mti', mti, 'avm', avm, 'taus', taus, 'rtests', rtests, 'errBCavm', zeros(Tl*kl,2+mti),...
        'errTrueF', zeros(Tl*kl,2+mti), 'nnz', zeros(Tl*kl,2), 'perc', zeros(Tl*kl,2), 'errSol', zeros(Tl*kl,2+mti), 'errBCcol', ...
        zeros(Tl*kl,2+mti), 'compresErr', zeros(Tl*kl,2), 'timeSol', zeros(Tl*kl,2+mti), 'nbIter', zeros(Tl*kl,mti), 'timeA', ...
        zeros(Tl*kl,2), 'ks', ks, 'field', zeros(70), 'errInt', zeros(Tl*kl,2+mti) );
else
    v = struct('perc', zeros(Tl*kl,2), 'errSol', zeros(Tl*kl,2+mti), 'mti', mti, 'avm', avm, 'taus', taus, 'rtests', rtests, 'errBCavm', ...
        zeros(Tl*kl,2+mti), 'timeSol', zeros(Tl*kl,2+mti), 'errBCcol', zeros(Tl*kl,2+mti), 'compresErr', zeros(Tl*kl,2), ...
        'timeA', zeros(Tl*kl,2), 'ks', ks, 'errInt', zeros(Tl*kl,2+mti), 'conds', zeros(Tl*kl,2), 'nbIter', zeros(Tl*kl,mti) );
end

%% Computations
for ki = 1:kl
    par.k = ks(ki); par.N = par.ppw*par.k;
    par.t = linspace(0,1,par.N + 1); % The knots of the periodic spline;
    par.colltau = par.t(1:par.N);
    
    A1 = zeros(par.N);
    counter = 0;
    tic;
    prevToc = toc;
    for i = 1:par.N
        if (toc-prevToc > 5)
            prevToc = toc;
            display(['k = ', num2str(par.k), ', ' num2str( (i-1)/par.N,'%7.3f'), ' = part of matrix, estim sec left = ' num2str(toc*(par.N-i+1)/(i-1)) ])
        end
        % 	parfor i=1:par.N % Possibly use this instead of a sequential loop
        A1(i,:) = collRowQBF(i,par);
    end
    ix = Tl*(ki-1)+1:Tl*ki;
    v.timeA(ix,1) = toc;
    
    for ti = 1:Tl
        T = Ts(ti);
        A2 = zeros(par.N);
        tic;
        prevToc = toc;
        for i = 1:par.N
            if (toc-prevToc > 5)
                prevToc = toc;
                display(['k = ', num2str(par.k), ', ' num2str( (i-1)/par.N,'%7.3f'), ' = part of compr. matrix, estim sec left = ' ...
                    num2str(toc*(par.N-i+1)/(i-1)) ])
            end
            tc = par.colltau(i);
            tbounds = [tc-T; tc+T; T*percDecay; T*percDecay];
            % Let windRow fix bounds as well as window plateau
            if     tc+T < 0.5-tc+factShad*T
                if tc-T > 0.5-tc-factShad*T
                    tbounds = [tbounds, [tc-T; 0.5-tc+factShad*T; T*percDecay; factShad*T*percDecay]];
                else
                    tbounds = [tbounds, [0.5-tc-factShad*T; 0.5-tc+factShad*T; factShad*T*percDecay; factShad*T*percDecay]];
                end
            end % No elseif because both are correct for large T
            if tc-T > 1.5-tc-factShad*T
                if tc+T < 1.5-tc+factShad*T
                    tbounds = [[1.5-tc-factShad*T; tc+T; factShad*T*percDecay; T*percDecay], tbounds];
                else
                    tbounds = [[1.5-tc-factShad*T; 1.5-tc+factShad*T; factShad*T*percDecay; factShad*T*percDecay], tbounds];
                end
            end
            tbounds(1,:) = tbounds(1,:) - eps; % Make sure that tbounds(1,:) < tbounds(2,:).
            tbounds(2,:) = tbounds(2,:) + eps;
            warning('off', 'windRow:removeGreen'); % When T is too small, the matrix will loose the Green singularity and become singular.
            A2(i,:) = windRow(i,par,tbounds);
        end
        ix = Tl*(ki-1)+ti;
        v.timeA(ix,2) = toc;
        v = validate(A1,A2,par,v,ix);
    end
end

%% Plots
if varykT
    clearvars A1 A2
    save fixedWindows.mat
    % Afterwards, do:
    if 0
        clearvars;
        load fixedWindows.mat;
        
        fs = 27;
        l = {'*', '+', 'o', '.', 'x', 's', 'd', '>', '<', 'p', 'h', 'v', '^'};
        c = 'bgrcmk';
        fss = 'Fontsize'; fs = 22;
        lws = 'LineWidth'; lw = 5;
        labs = [repmat('k = ', length(ks),1), num2str(ks')];
        matr = reshape(v.errSol(:,2),Tl,kl);
        matr(abs(matr) > 1) = 1;
        figure;
        for ki = 1:kl
            semilogy(Ts, matr(:,ki), [l{ki} c(mod(ki-1,length(c))+1)], 'MarkerSize',10); hold on;
        end
        legend(labs);
        xlabel('T','FontSize',fs);
        ylabel('$|| \tilde{c}-$c$ ||/|| $c$ ||$','interpreter','latex',fss,fs);
        set(gca,fss,fs);
    end
else
    plotVal(struct(),A2);
    % Test using a block window.
    A3 = A1;
    A3(A2==0) = 0;
    v3 = validate(A1, A3, par, struct('perc', zeros(Tl*kl,2), 'errSol', zeros(Tl*kl,2+mti), 'mti', mti, 'avm', avm, 'taus', taus, ...
        'rtests', rtests, 'errBCavm', zeros(Tl*kl,2+mti), 'timeSol', zeros(Tl*kl,2+mti), 'errBCcol', ...
        zeros(Tl*kl,2+mti), 'compresErr', zeros(Tl*kl,2),  'ks', ks, 'errInt', zeros(Tl*kl,2+mti), 'conds', zeros(Tl*kl,2) ), ix)
    % Check the sizes of some contributions:
    b = par.bc(par.k,par.par(par.colltau));
    c1 = A1\b;
    for row = [round(par.N/8), round(par.N*5/9)] % Some row in the shadow region and in the illuminated region
        gs = par.colltau(row);
        greenSingCtrb = (chi(par.colltau, gs-T/4, gs-T/8, gs+T/8, gs+T/4,1).*A1(row,:) )*c1;
        figure;
        plot(par.colltau, [chi(par.colltau, gs-T/4, gs-T/8, gs+T/8, gs+T/4,1); abs(A3(row,:).*transpose(c1))/max(abs(A3(row,:).*transpose(c1))) ;...
            real(A3(row,:).*transpose(c1))/max(abs(A3(row,:).*transpose(c1)));...
            chi(par.colltau, gs-5*T/4, gs-9*T/8, gs-7*T/8, gs-T*3/4,1); chi(par.colltau, 0.5-gs-0.2, 0.5-gs-0.15, 0.5-gs+0.15, 0.5-gs+0.2,1)]);
        endPtCtrb = (chi(par.colltau, gs-5*T/4, gs-9*T/8, gs-7*T/8, gs-T*3/4,1).*A3(row,:) )*c1;
        statPtCtrb = (chi(par.colltau, 0.5-gs-0.2, 0.5-gs-0.15, 0.5-gs+0.15, 0.5-gs+0.2,1).*A1(row,:) )*c1;
        display(['For row ' num2str(row) ', the absolute value of the contribution of the green singularity is ' ...
            num2str(abs(greenSingCtrb)) ', of the SP is ' num2str(abs(statPtCtrb)) ' and of the endpoint is ' num2str(abs(endPtCtrb))]);
        display(['The sum of the former is ' num2str(greenSingCtrb + statPtCtrb) ' while b(' num2str(row) ') is ' num2str(b(row)) ]);
    end
end
v
