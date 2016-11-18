% Compute the solution of multiple scattering obstacles using a known phase.
if 0
%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

ks = 2^7;
obsts = 5;
bc = 1;
% bcsh = pi/2; % shift in boundary condition

printtoc = 3;
kl = length(ks);
nbOb = length(obsts);

avm = 100; % Number of random taus to average BC over
% v = struct('mti', mti, 'avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(nbOb*kl,2+mti),...
%     'perc', zeros(nbOb*kl,2), 'errSol', zeros(nbOb*kl,2+mti),  'errInt', zeros(nbOb*kl,2+mti), ...
%     'timeSol', zeros(nbOb*kl,2+mti), 'nbIter', zeros(nbOb*kl,mti), 'timeA', zeros(nbOb*kl,4), 'ks', ks, 'field', zeros(70) );
v = struct('avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(nbOb*kl,2), 'errSol', zeros(nbOb*kl,2),  ...
    'errInt', zeros(nbOb*kl,2), 'timeA', zeros(nbOb*kl,4), 'ks', ks); %, 'field', zeros(70) );

%% Computations
for oi = 1:length(obsts)
    obstacle = obsts(oi);
    start = now;
    
    for ki = 1:kl
        idx = (oi-1)*kl+ki;
        par = getObst(obstacle); % Reset par
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
        
        if bc == 4
            par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/2-bcsh);sin(-pi/2-bcsh)]);
        elseif bc == 3
            par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/2);sin(-pi/2)]);
        elseif bc == 2
            par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/4);sin(-pi/4)]);
        elseif bc == 1
            par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(0);sin(0)]);
        end
        
        if isfield(par,'phase')
            par = rmfield(par, 'phase'); % no phase when computing A1
            %             par = rmfield(par, 'quadR'); % don't need additional points when computing A1
        end
        overs1 = 1;
        
        %% Computating full solution
%         A1 = zeros(par.N); tic;
%         A1 = zeros(overs1*par.N); 
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
%                     A1(i,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N);
                else
                    for os = 1:overs1
                        A1(overs1*i-overs1+os,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxos(:,os) );
                    end
%                     A1(i,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N,  [], collx);
                end
            end
        end
        v.timeA(idx,1) = toc;
        b = zeros(par.N,1);
        for obst = 1:length(par.obsts)
            b(par.r(1,obst):par.r(2,obst)) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).colltau));
        end
%         c1 = A1\b; % The solution is needed for computing the correlations.
        c1 = A1\b1; % The solution is needed for computing the correlations.
        return
        
        % [R, sigma,obbounds] = calcCorr(par, c1, Tcor, percDecay); 
        % [R, sigma,obbounds] = calcCorr(par, c1, 0.1, 1, [10, now], A1); 
        
        %% Compute A2
%         fN = 48;
        fN = round(par.N/length(par.obsts));
%         fN = round(1280/32);
%         fN = round(1280/8);
        ft = linspace(0,1,fN+1); % The knots of the periodic spline;
        fco = ft(1:fN);
%         par.quadR = 25;
%         par.quadR = 5;
        par.quadR = max(2, 6*round(par.N/fN) );
%         for obst = 1:length(par.obsts)
%             par.obsts(obst).quadR = max(2, 6*round(par.N/fN) );
%         end
        if obstacle == 5 && bc == 1
            nrPh = 2; % Number of phases
%             overs = 2;
            overs = 7;
        elseif obstacle == 5 && bc == 3
            nrPh = 1;
        elseif obstacle == 5 && bc == 4
            nrPh = 1;
        end
%         A2 = zeros(length(par.obsts)*nrPh*fN);
        A2 = zeros(length(par.obsts)*overs*fN,length(par.obsts)*nrPh*fN);
        
        tmp = par;
        par.N = 0;
        par.r = zeros(2,length(par.obsts)); % ranges
        for obst = 1:length(par.obsts)
            par.obsts(obst).N = fN; %par.obsts(obst).ppw*fN;
            par.r(1,obst) = par.N+1;
            par.N = par.N + par.obsts(obst).N;
            par.r(2,obst) = par.N;
            
            par.obsts(obst).t = ft; % The knots of the periodic spline;
            par.obsts(obst).colltau = fco;
            
%             par.obsts(obst).quadR = max(2, 6*round(par.N/fN) );
            par.obsts(obst).quadR = par.quadR;
            par.obsts(obst).fco = fco;
            par.obsts(obst).fN = fN;
            if obstacle == 5 && bc == 1
                par.obsts(obst).phase = @(t) transpose(transpose(par.obsts(obst).par(t) )*[cos(0);sin(0)]);
            elseif obstacle == 5 && bc == 4
                par.obsts(obst).phase = @(t) transpose(transpose(par.obsts(obst).par(t) )*[cos(-pi/2-bcsh);sin(-pi/2-bcsh)]);
%                 par.obsts(obst).phase = @(t) transpose(transpose(par.obsts(obst).par(t - 0.5*(t>0.5) ) )*[cos(-pi/2-bcsh);sin(-pi/2-bcsh)]); %Still oscillations
            elseif obstacle == 5 && bc == 3
                par.obsts(obst).phase = @(t) transpose(transpose(par.obsts(obst).par(t) )*[cos(-pi/2);sin(-pi/2)]); % Or different in shadow
%                 par.obsts(obst).phase = @(t) transpose(transpose(par.obsts(obst).par(t - 0.5*(t>0.5)) )*[cos(-pi/2);sin(-pi/2)]); % but
%                 gives oscillations in shadow probably due to nonanalyticity
            end
            par.obsts(obst).hs = (par.obsts(obst).colltau(2) -par.obsts(obst).colltau(1) )/2;
        end
        tic;
        prevToc = toc;
        nN = length(par.obsts)*fN;
        
        b2 = nan*zeros(size(A2,1),1);
        for i = 1:nN
            if (toc-prevToc > printtoc)
                prevToc = toc;
                display([num2str( (i-1)/nrPh/nN,'%7.3f'), '=A2%, first, est. # sec. left for A2=' num2str(toc*(nrPh*nN-i+1)/(i-1)) ])
            end
            obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
            collx = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) );
            collxsh4 = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) + par.obsts(obsin).hs/2);
            collxsh = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) + par.obsts(obsin).hs);
            collxsh34 = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) + par.obsts(obsin).hs*3/2);
            collxos = zeros(2,overs);
            for os = 1:overs
                collxos(:,os) = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) + par.obsts(obsin).hs*2/overs*(os-1));
            end
            for obst = 1:length(par.obsts)
                if obst == obsin && nrPh == 1
                    A2(i,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N);
                    b2(i) = par.bc(par.k, collx);
%                 elseif obst == obsin && nrPh == 2
                elseif obst == obsin
                    for os = 1:overs
                        A2(overs*i-overs+os,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N, [], [], ...
                            par.obsts(obst).hs*2/overs*(os-1) );
                        b2(overs*i -overs+os) = par.bc(par.k, collxos(:,os));
                    end
%                 elseif obst == obsin && overs == 2
%                     A2(2*i-1,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N);
% %                     A2(2*i,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N);
%                     A2(2*i,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N, [], [], par.obsts(obst).hs);
%                     b2(2*i - 1) = par.bc(par.k, collx);
%                     b2(2*i) = par.bc(par.k, collxsh);
%                 elseif obst == obsin && overs == 4
%                     A2(4*i-3,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N);
%                     A2(4*i-2,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N, [], [], par.obsts(obst).hs/2);
%                     A2(4*i-1,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N, [], [], par.obsts(obst).hs);
%                     A2(4*i,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N, [], [], par.obsts(obst).hs*3/2);
%                     b2(4*i - 3) = par.bc(par.k, collx);
%                     b2(4*i - 2) = par.bc(par.k, collxsh4);
%                     b2(4*i - 1) = par.bc(par.k, collxsh);
%                     b2(4*i) = par.bc(par.k, collxsh34);
                elseif nrPh == 1
                    A2(i,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collx);
                else
                    for os = 1:overs
                        A2(overs*i-overs+os,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxos(:,os) );
                    end
%                 elseif overs == 2
%                     A2(2*i-1,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collx);
%                     A2(2*i,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxsh);
% %                     A2(nN + i,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxsh);
%                 elseif overs == 4
%                     A2(4*i-3,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collx);
%                     A2(4*i-2,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxsh4);
%                     A2(4*i-1,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxsh);
%                     A2(4*i,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxsh34);
%                 else
%                     error('Check overs.');
                end
            end
        end
        
        if nrPh > 1
%             par.obsts(1).phase = @(t) transpose(transpose(par.obsts(1).par(t) )*[cos(0);sin(0)]);
%             par.obsts(2).phase = @(t) transpose(transpose(par.obsts(2).par(t) )*[cos(0);sin(0)]);
            par.obsts(1).phase = @(t) sqrt( sum( (par.obsts(1).par(t) - repmat(par.obsts(2).par(3/4), 1, length(t) ) ).^2, 1) );
            par.obsts(2).phase = @(t) sqrt( sum( (par.obsts(2).par(t) - repmat(par.obsts(1).par(1/4), 1, length(t) ) ).^2, 1) );
%             par.obsts(1).phase = @(t) sqrt( -sum( (par.obsts(1).par(t) - repmat(par.obsts(2).par(3/4), 1, length(t) ) ).^2, 1) );
%             par.obsts(1).phase = @(t) -sqrt( sum( (par.obsts(1).par(t) - repmat(par.obsts(2).par(3/4), 1, length(t) ) ).^2, 1) );
%             par.obsts(2).phase = @(t) -sqrt( sum( (par.obsts(2).par(t) - repmat(par.obsts(1).par(1/4), 1, length(t) ) ).^2, 1) );
%             par.obsts(1).phase = @(t) arrayfun(@(tt) norm(par.obsts(1).par(tt) - par.obsts(2).par(3/4) ), t);
%             par.obsts(2).phase = @(t) arrayfun(@(tt) norm(par.obsts(2).par(tt) - par.obsts(1).par(1/4) ), t);
%             par.obsts(2).phase = @(t) transpose(transpose(par.obsts(obst).par(t - 1/4) )*[cos(0);sin(0)]);
            %         par.obsts(1).phase = @(t) transpose(transpose(par.obsts(obst).par(t - 1/4) )*[cos(0);sin(0)]);
            %         par.obsts(2).phase = @(t) transpose(transpose(par.obsts(obst).par(t + 1/4) )*[cos(0);sin(0)]);
            
            for i = 1:nN
                if (toc-prevToc > printtoc)
                    prevToc = toc;
                    display([num2str( 0.5+(i-1)/2/nN,'%7.3f'), '=A2%, second, est. # sec. left for A2=' num2str(toc*(nN-i+1)/(i-1+nN)) ])
                end
                obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
                collx = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) );
                collxsh = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) + par.obsts(obsin).hs);
                for obst = 1:length(par.obsts)
                    if obst == obsin
                        for os = 1:overs
                            A2(overs*i-overs+os, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N, [], [], ...
                                par.obsts(obst).hs*2/overs*(os-1) );
                        end
                    else
                        for os = 1:overs
                            A2(overs*i-overs+os, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxos(:,os) );
                        end
                    end
%                     if (obst == obsin) && overs == 2
%                         A2(2*i-1, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N);
%                         A2(2*i, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N, [], [], par.obsts(obst).hs);
%                     elseif (obst == obsin) && (overs == 4)
%                         A2(4*i-3, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N);
%                         A2(4*i-2, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N, [], [], par.obsts(obst).hs/2);
%                         A2(4*i-1, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N, [], [], par.obsts(obst).hs);
%                         A2(4*i, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N, [], [], par.obsts(obst).hs*3/2);
%                     elseif overs == 2
%                         A2(2*i-1, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collx);
%                         A2(2*i, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxsh);
%                         %                     A2(nN + i,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxsh);
%                     elseif overs == 4
%                         A2(4*i-3, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collx);
%                         A2(4*i-2, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxsh4);
%                         A2(4*i-1, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxsh);
%                         A2(4*i, nN + (par.r(1,obst):par.r(2,obst))) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxsh34);
%                     else
%                         error('Check overs');
%                     end
                end
            end
%             
%             b2t = nan*zeros(size(A2,1),1);
%             for obst = 1:length(par.obsts)
%                 b2t(2*(par.r(1,obst):par.r(2,obst))-1) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).colltau));
%                 %             b(nN + par.r(1,obst):par.r(2,obst)) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).colltau + par.obsts(obst).hs));
%                 b2t(2*(par.r(1,obst):par.r(2,obst))) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).colltau + par.obsts(obst).hs));
%             end
        else
%             b2 = b;
        end
%         b2 = par.bc(par.k,par.par(par.fco));
        c2 = A2\b2;
%         erSol = norm(interp1([par.fco (par.fco(1)+1)], [c2; c2(1)], par.colltau').*transpose(exp(1i*par.k*par.phase(par.colltau) ))-c1)/norm(c1)
        par.N = tmp.N;
        par.fr = par.r;
        par.r = tmp.r;
        for obst = 1:length(par.obsts)
            par.obsts(obst).N = tmp.obsts(obst).N;
            par.obsts(obst).t = tmp.obsts(obst).t; 
            par.obsts(obst).colltau = tmp.obsts(obst).colltau;
        end
        
        %% Testing
%         ct = ones(size(c2));
        ct = zeros(size(c2));
%         ct(1) = 1;
%         ct(end) = 1;
        col = par.obsts(1).colltau;
%         ct(round(length(c2)/8):3*round(length(c2)/8)) = cos(2*pi*col(round(length(col)/4):3*round(length(col)/4)) );
%         ct(round(length(c2)/16):round(3*length(c2)/16)) = cos(2*pi*col(round(length(col)/4):3*round(length(col)/4)) );
%         shft = 1/2; % 1/4;
        shft = 3/4-1/16;
%         ct(round(length(c2)*(shft+1/16)):round(length(c2)*(shft+3/16) )) = cos(2*pi*col(round(length(col)/4):3*round(length(col)/4)) );
%         bt = A2*ct;
%         figure; plot([real(bt) imag(bt)])
%         bco = par.bc(par.k, par.obsts(1).par(col));
%         figure; plot(col, bco);
        
        %% finishing
%         figure; plot(par.colltau, par.phase(par.colltau)); title('Phase')
        figure; spectrogram(c1, round(length(c1)/16),[],[], 1/(par.obsts(1).colltau(2)-par.obsts(1).colltau(1)), 'centered'); title(['Spectrogram ' ])%num2str(bcsh)])
        figure; plot( [real(c1) imag(c1)]); legend('Re(c1)','Im(c1)')
        figure; plot( [real(c2) imag(c2)]); legend('Re(c2) using fco','Im(c2) using fco')
        v = validate(A1,A2,par,v,idx)
%         plotVal(v)
        return

        %% angle of part where only phase indep of inc wave
        collsignal = par.obsts(2).colltau((length(par.obsts(2).colltau)*0.8-1):end);
        signal = c1(end-round(length(c1)*0.1+1):end);
%         figure; plot(collsignal, [real(signal) imag(signal)]);
        l = length(signal);
        ph = nan*signal;
        for i = l:-1:1
            ph(i) = angle(signal(i));
            if (i ~= l) && (abs(ph(i) - ph(i+1)) > 5)
                ph(i) = ph(i) -2*pi*round( (ph(i)-ph(i+1))/2/pi);
%             elseif (i ~= l) && (ph(i) - ph(i+1) < 5)
%                 ph(i) = ph(i) +2*pi*round( (ph(i)-ph(i+1))/2/pi);
            end
%             elseif phase(i) - phase(is(nr-1)) > 5
%         phase(i) = phase(i) -2*pi*round( (phase(i)-phase(is(nr-1)))/2/pi);
%     elseif phase(is(nr-1)) - phase(i) > 5
%         phase(i) = phase(i) +2*pi*round( (phase(is(nr-1))-phase(i))/2/pi);
        end
        figure; plot(collsignal, ph/par.k)

        %% Selfmade spectrogram without windows to get ifft
%         signal = c1(end-length(c1)/2+1:end-length(c1)/4);
%         collsignal = par.obsts(2).colltau((length(par.obsts(2).colltau)/2+1):end);
%         signal = c1(end-length(c1)/4+1:end);
        collsignal = par.obsts(2).colltau;
        signal = c1(end-length(c1)/2+1:end);
%         figure; plot(collsignal, [real(signal) imag(signal)]);

        nrTim = 16;
        nrFreq = ceil(length(signal)/nrTim);
        S = zeros(nrFreq, nrTim);
        for nt = 1:nrTim
            S(:,nt) = fft(signal( (nrFreq*(nt-1)+1):(nrFreq*nt) ) );
        end
        for nt = 1:0%2
            figure;
            pcolor(abs(S) );
            hold on
            shading interp;
            xlabel('time');
            ylabel('frequency');
            set(gca,'FontSize',20);
            colorbar();
            if bc == 3
                S(1:10, 1:5) = 0;
                S(end-9:end, end-4:end) = 0;
            elseif bc == 1;
                S(1,1:10) = 0;
                S(2:3,1:11) = 0;
                S(4,1:12) = 0;
                S(5,1:13) = 0;
                S(6,1:14) = 0;
                S(7:74,:) = 0;
                S(75,1:6) = 0;
                S(76:77, 1:7) = 0;
                S(78, 1:8) = 0;
                S(79:80, 1:9) = 0;
            end
        end
        reconstr = nan*signal;
        for nt = 1:nrTim
            reconstr( (nrFreq*(nt-1)+1):(nrFreq*nt) ) = ifft(S(:,nt) );
        end
%         figure; 
%         plot([real(signal) imag(signal) real(reconstr) imag(reconstr)]);
%         legend('Re signal', 'Im signal', 'Re reconstr', 'Im reconstr');
%         figure; plot([real(signal) real(reconstr)]); legend('Re signal','Re reconstr');
%         figure; plot([imag(signal) imag(reconstr)]); legend('Im signal', 'Im reconstr');

        par.obsts(2).phase = @(t) sqrt( sum( (par.obsts(2).par(t) - repmat(par.obsts(1).par(1/4), 1, length(t) ) ).^2, 1) );
        deconv = reconstr.*exp(-1i*par.k*transpose(par.obsts(2).phase(collsignal)));
%         figure; plot(collsignal, [real(deconv) imag(deconv)]);
%         figure; spectrogram(deconv, round(length(deconv)/16),[],[], collsignal(2)-collsignal(1), 'centered'); 
    end
%     v = validate(A1,A2,par,v,idx)
    display([num2str(oi) ' = oi, ki = ' num2str(ki) ', now is ' datestr(now) ', expected end ' datestr(start + ...
        (now-start)*sum(ks.^2)/sum(ks(1:ki).^2)*( length(obsts)-oi + 1) )  ]);
end

% return

%% Construct c3 from presumed c2
c3 = 0*c1;
par.obsts(1).phase
for obst = 1:length(par.obsts)
%     c3(par.r(1,obst):par.r(2,obst)) = transpose(exp( (par.obsts(obst).colltau - (2*obst-1)/4).^2 ));
    c3(par.r(1,obst):par.r(2,obst)) = 80*transpose(exp( -abs(sin(2*pi*par.obsts(obst).colltau) - sin(pi*(2*obst-1)/2) )*4.5 ).*...
        exp(1i*par.k*par.obsts(obst).phase(par.obsts(obst).colltau) ));
end

for obst = 1:length(par.obsts)
    c3(par.r(1,obst):par.r(2,obst)) = c3(par.r(1,obst):par.r(2,obst)) + 280*transpose(exp( -abs(cos(2*pi*par.obsts(obst).colltau) +1)*3.1 ).*...
        exp(1i*par.k*[1, 0]*par.obsts(obst).par(par.obsts(obst).colltau) ));
end

figure; plot([real(c3) imag(c3)]);

%% Construct presumed c2
c4 = 0*c2;
% par.obsts(1).phase

for obst = 1:length(par.obsts)
    c4(par.r(1,obst):par.r(2,obst)) = 1i*280*transpose(exp( -abs(cos(2*pi*par.obsts(obst).colltau) +1)*3.1 ));
end

for obst = 1:length(par.obsts)
    c4(par.r(2,length(par.obsts)) + (par.r(1,obst):par.r(2,obst))) = 80*transpose(exp( -abs(sin(2*pi*par.obsts(obst).colltau) - sin(pi*(2*obst-1)/2) )*4.5) );
end
% figure; plot([real(c4) imag(c4)]); legend('Re c4', 'Im c4')
b4 = A2*c4;
figure; plot([real(b4) imag(b4)]); legend('Re b4', 'Im b4')

%% Remove columns associated to c2-entries which are expected to be zero
A3 = A2;
% Start from the end because the indices will change
A3(:, (size(A2,2)*3/4 + 1):(size(A2,2)*7/8)) = [];
A3(:, (size(A2,2)*1/2 + 1):(size(A2,2)*5/8)) = [];
A3(:, (size(A2,2)/4 + find(par.obsts(2).colltau > 0.85, 1, 'first')):size(A2,2)/2) = [];
A3(:, (size(A2,2)/4 + 1):(size(A2,2)/4 + find(par.obsts(2).colltau > 0.15, 1, 'first'))) = [];
A3(:, find(par.obsts(1).colltau > 0.85, 1, 'first'):size(A2,2)/4) = [];
A3(:, 1:find(par.obsts(1).colltau > 0.15, 1, 'first')) = [];

c5 = A3\b2;

%% Add rows associated to c2-entries which are expected to be zero

% Start from the end because the indices will change
idxs = [(size(A2,2)*3/4 + 1):(size(A2,2)*7/8), (size(A2,2)*1/2 + 1):(size(A2,2)*5/8), ...
    (size(A2,2)/4 + find(par.obsts(2).colltau > 0.85, 1, 'first')):size(A2,2)/2, (size(A2,2)/4 + 1):(size(A2,2)/4 + find(par.obsts(2).colltau > 0.15, 1, 'first') ), ...
    find(par.obsts(1).colltau > 0.85, 1, 'first'):size(A2,2)/4, 1:find(par.obsts(1).colltau > 0.15, 1, 'first')];
A4 = [A2; zeros(length(idxs), size(A2,2))];
for nix = 1:length(idxs)
    A4(size(A2,1)+nix, idxs(nix)) = 1;
end
b6 = [b2; zeros(length(idxs),1)];
c6 = A4\b6;

%% remove FFT components
cut = ceil(length(c1)/par.obsts(1).ppw/length(par.obsts)/4);
% cut = 8;
tmp = fft(c1.*transpose(exp(-1i*par.k*[1, 0]*[par.obsts(1).par(par.obsts(1).colltau), par.obsts(2).par(par.obsts(2).colltau)] ) ));
tmp((cut+1):(end-cut)) = 0;
c1inc = ifft(tmp).*transpose(exp(1i*par.k*[1, 0]*[par.obsts(1).par(par.obsts(1).colltau), par.obsts(2).par(par.obsts(2).colltau)] ) );
figure; plot([real(c1inc) imag(c1inc)]); title('c1inc')
resi = c1-c1inc;
figure; plot([real(resi) imag(resi)]); title('residual')
figure; spectrogram(resi, round(length(resi)/16),[],[], 1/(par.obsts(1).colltau(2)-par.obsts(1).colltau(1)), 'centered');

phas = exp(1i*par.k*[sqrt( sum( (par.obsts(1).par(par.obsts(1).colltau) - repmat(par.obsts(2).par(3/4), 1, length(par.obsts(1).colltau) ) ).^2, 1) ), ...
    sqrt( sum( (par.obsts(2).par(par.obsts(2).colltau) - repmat(par.obsts(1).par(1/4), 1, length(par.obsts(2).colltau) ) ).^2, 1) ) ]);
tmp = resi.*transpose(1./phas);
% tmp = c1.*transpose(1./phas);
tmpfft = fft(tmp);
tmpfft((cut+1):(end-cut)) = 0;
c1sp = ifft(tmpfft).*transpose(phas);
figure; plot([real(c1sp) imag(c1sp)]); title('c1sp');
figure; plot([real(ifft(tmpfft)) imag(ifft(tmpfft))]); title('ifft(tmpfft)');
figure; spectrogram(c1sp, round(length(c1sp)/16),[],[], 1/(par.obsts(1).colltau(2)-par.obsts(1).colltau(1)), 'centered');

%% FFT of slices
slice = 180:1100;
cut = 16;
c1slill = c1(slice);
phas = transpose(exp(1i*par.k*[1, 0]*par.obsts(1).par(par.obsts(1).colltau(slice) )) );
ffttmp = fft(c1slill./phas);
ffttmp((cut+1):(end-cut)) = 0;
c1inc = ifft(ffttmp).*phas;
figure; plot([abs(ffttmp) abs(fft(c1slill./phas))]) 
figure; plot([real(ifft(ffttmp)) imag(ifft(ffttmp))]); title('c1incLF')
figure; plot([real(c1inc) imag(c1inc)]); title('c1inc')


%% Substract solution from one circle with phase shift
load c1OneCirc.mat
slice = 1:length(par.obsts(1).colltau); 
% figure; plot(real(c1(slice)./c1OneCirc));
c1s = c1(slice);
figure; plot([real(c1s) imag(c1s)]); title('c1slice');
rat = mean(c1s(600:1150)./c1OneCirc(600:1150))
c1sp = c1s - c1OneCirc*rat;
% c1iLF = c1s.*transpose(exp(-1i*par.k*[1, 0]*par.obsts(1).par(par.obsts(1).colltau)));
c1iLF = c1OneCirc*rat.*transpose(exp(-1i*par.k*[1, 0]*par.obsts(1).par(par.obsts(1).colltau)));
figure; plot([real(c1iLF) imag(c1iLF)]); title('c1illf');
figure; plot([real(c1sp) imag(c1sp)]); title('c1sp');
figure; spectrogram(c1sp, round(length(c1sp)/16),[],[], 1/(par.obsts(1).colltau(2)-par.obsts(1).colltau(1)), 'centered');

%% reaosifdjh
phas = exp(1i*par.k*sqrt( sum( (par.obsts(1).par(par.obsts(1).colltau) - repmat(par.obsts(2).par(3/4), 1, length(par.obsts(1).colltau) ) ).^2, 1) ));
% phas = exp(1i*par.k*sqrt( sum( (par.obsts(1).par(par.obsts(1).colltau) - repmat(par.obsts(2).par(3/4+0.15), 1, length(par.obsts(1).colltau) ) ).^2, 1) ));
tmp = c1sp.*transpose(1./phas);
figure; plot([real(tmp) imag(tmp)]); title('tmp');
figure; plot(abs(fft(tmp))); title('ffttmp')

c1sp2 = c1(slice+length(par.obsts(1).colltau)) - c1OneCirc*rat;
figure; plot([real(c1sp2) imag(c1sp2)]); title('c1sp2')
tmp2 = c1sp2.*transpose(exp(-1i*par.k*sqrt( sum( (par.obsts(2).par(par.obsts(2).colltau) - repmat(par.obsts(1).par(1/4), 1, length(par.obsts(2).colltau) ) ).^2, 1) )));

%% Check using the lowfreq part
% c2t = nan*c2;
% c2t = [c1OneCirc*rat; c1OneCirc*rat; tmp; flipud(tmp)];
% c2t = [c1iLF; c1iLF; tmp; flipud(tmp)];
c2t = [c1iLF; c1iLF; tmp; tmp2];
b2t = A2*c2t;
figure; plot([real(b2t) imag(b2t)]); title('b2t')
figure; plot([real(c2t) imag(c2t)]); title('c2t')
errb = norm(b2t-b2)/norm(b2)

%% Coincident phases?
phsp = sqrt( sum( (par.obsts(1).par(par.obsts(1).colltau) - repmat(par.obsts(2).par(3/4), 1, length(par.obsts(1).colltau) ) ).^2, 1) );
phil = [1, 0]*par.obsts(1).par(par.obsts(1).colltau);
figure; plot(par.obsts(1).colltau, [phsp; phil])
figure; plot(par.obsts(1).colltau, cos(par.k*[phsp; phil])); title('cos(k*ph)');

%% eigenvalue decomp
[V, eigA2, W] = eig(A2);

%% plot
figure; plot(log(abs(diag(eigA2))))
figure; plot([real(V(:,end)) imag(V(:,end))])
figure; plot([real(W(:,1)) imag(W(:,1))])

%% Remove components
rem = 100;
% rem = 4700;
% c2f = c2;
b2f = b2;
for ri = 1:rem
%     c2f = c2f - V(:,end-ri+1)*(transpose(V(:,end-ri+1))*c2);
%     c2f = c2f - V(:,end-ri+1)*(transpose(V(:,end-ri+1))*c2f);
    b2f = b2f - V(:,end-ri+1)*(transpose(V(:,end-ri+1))*b2);
end
% eigenvectors normal but not orthogonal wrt each other?
% [norm(c2), norm(c2f), norm(c2t), rem]
display([num2str(rem) ' components are removed'])
c2f = A2\b2f;
figure; plot([real(c2f) imag(c2f)]); title('c2f');


%% Ray tracing by partitioned block 
A11 = A1(1:size(A1,1)/2, 1:size(A1,2)/2);
A21 = A1(1:size(A1,1)/2, (1+size(A1,2)/2):end); % Action of density of obst 2 on location of obst 1
A12 = A1((size(A1,1)/2+1):end, 1:size(A1,2)/2);
A22 = A1((size(A1,1)/2+1):end, (1+size(A1,2)/2):end);

c1true1 = c1(1:length(c1)/2);
c1true2 = c1(length(c1)/2+1:end);

nbrefl = 8;
bsl1 = zeros(size(A1,1)/2, nbrefl);
bsl2 = zeros(size(A1,1)/2, nbrefl);
csl1 = zeros(size(A1,2)/2, nbrefl);
csl2 = zeros(size(A1,2)/2, nbrefl);

bsl1(:,1) = b(1:size(A1,1)/2);
bsl2(:,1) = b((size(A1,1)/2+1):end);
csl1(:,1) = A11\bsl1(:,1);
csl2(:,1) = A22\bsl2(:,1);

resid = cell(2,nbrefl);
resid{1,1} = c1true1 - sum(csl1,2);
figure;
plot(par.obsts(1).colltau, real(csl1(:,1))/norm(csl1(:,1)) ); title('csl1 refl normal')
% plot(par.obsts(1).colltau, real(resid{1,1})); title('resid')
% plot(par.obsts(1).colltau, real(sum(csl1,2))); title('sum')
hold on

for refl = 2:nbrefl
%     for refl = 1:size(csl1,2);
%     if refl == 1
%         bsl1(:,1) = b(1:size(A1,1)/2);
%         bsl2(:,1) = b((size(A1,1)/2+1):end);
%     else
%     bsl1(:,refl) = A21*csl2(:,refl-1);
%     bsl2(:,refl) = A12*csl1(:,refl-1);
    bsl1(:,refl) = -A21*csl2(:,refl-1);
    bsl2(:,refl) = -A12*csl1(:,refl-1);
    csl1(:,refl) = A11\bsl1(:,refl); % Or coupling matrices?
    csl2(:,refl) = A22\bsl2(:,refl);
    
    resid{1,refl} = c1true1 - sum(csl1,2);
    figure; plot(par.obsts(1).colltau, real(csl1(:,refl))/norm(csl1(:,refl)) );
%     plot(par.obsts(1).colltau, real(resid{1,refl}));
%     plot(par.obsts(1).colltau, real(sum(csl1,2)));
end
legend(num2str((1:nbrefl)'));
% figure; plot(real(csl1))
% resid = cell(2);
% resid{1} = c1true1 - sum(csl1,2);
% figure; plot(par.obsts(1).colltau, [real(resid{1}) imag(resid{1})]);
% hold on;
% resid{1} = c1true1 - csl1(:,1);
% plot(par.obsts(1).colltau, [real(resid{1}) imag(resid{1})]);

%% Eigenvalues big reflection matrix

bigm = A11\A21*(A22\A12);
[Vb, Db, Wb] = eig(bigm);
% [V1, D1, W1] = eig(A1);
figure; semilogy(abs(diag(Db)))
figure; plot(par.obsts(1).colltau, [real(Vb(:,1)) imag(Vb(:,1))])
figure; plot(par.obsts(1).colltau, [real(Wb(:,1)) imag(Wb(:,1))])

cnor = sqrt(sum(abs(csl1).^2, 1))
ratcnor = cnor(1:end-1)./cnor(2:end)

% load('c1OneCirc.mat', 'c1sin'); % see phasCircle.m

% c1slice = c1(1:length(c1)/2);
% evd = V(:,1:8)\c1slice
% evd = Vb(:,1:8)\c1sp

%% Project eigenvectors of bigm onto assumed phase
obst = 1;
phbigm = exp(1i*par.k*sqrt( sum( (par.obsts(obst).par(par.obsts(obst).colltau) - repmat(par.obsts(3-obst).par((2*obst-1)/4), 1, ...
    length(par.obsts(obst).colltau) ) ).^2, 1) ) );
lfbm = Vb(:,1)./transpose(phbigm);
figure; plot(par.obsts(1).colltau, [real(lfbm) imag(lfbm)])

end % if 0 at begin

%% project c1 onto osc basis fcts  with low fN
% fN = 55;
fN = 300;
nrPh = 3;
A = zeros(size(c1,1), nrPh*fN*length(par.obsts));

ft = linspace(0,1,fN+1); % The knots of the periodic spline;
fco = ft(1:fN);
hco = fco(2)-fco(1);

phi = @(tau, i) (tau > fco(max(1,i-1)) ).*(tau <= fco(i) ).*(tau-fco(max(1,i-1)))/hco + (tau > fco(i)).*(tau < fco(min(fN,i+1)) ).*(1+(fco(i)-tau)/hco) + ...
    (tau > fco(end) ).*(i == 1).*(tau-fco(end))/hco + (tau > fco(end) ).*(i == fN).*(1+(fco(end)-tau)/hco);
% test = par.obsts(1).colltau*nan;
i = fN;
% for j = 1:length(test)
%     test(j) = phi(par.obsts(1).colltau(j), i);
% end
test = phi(par.obsts(1).colltau, i);
% figure; plot(par.obsts(1).colltau, test);

for obst = 1:length(par.obsts)
%     for obi = 1:length(par.obsts)
        for i = 1:fN
%             A(:,(obst-1)*fN+i) = exp(1i*par.k*[1,0]*par.obsts(obst).par( );
%             A(length(par.obsts(obi).colltau)*(obi-1)+1:length(par.obsts(obi).colltau)*obi ,(obst-1)*fN+i) = ...
%                 exp(1i*par.k*[1,0]*par.obsts(obst).par(par.obsts(obi).colltau) );
%             A(length(par.obsts(obst).colltau)*(obst-1)+1:length(par.obsts(obst).colltau)*obst, (obst-1)*fN+i) = ...
%                 exp(1i*par.k*[1,0]*par.obsts(obst).par(par.obsts(obst).colltau) );
            A(length(par.obsts(obst).colltau)*(obst-1)+1:length(par.obsts(obst).colltau)*obst, (obst-1)*fN+i) = ...
                phi(par.obsts(obst).colltau, i).*exp(1i*par.k*[1,0]*par.obsts(obst).par(par.obsts(obst).colltau) );
        end
%     end
end
% now second phase:
for obst = 1:length(par.obsts)
    for i = 1:fN
        A(length(par.obsts(obst).colltau)*(obst-1)+1:length(par.obsts(obst).colltau)*obst, (length(par.obsts)+obst-1)*fN+i) = ...
            exp(1i*par.k*sqrt( sum( (par.obsts(obst).par(par.obsts(obst).colltau) - repmat(par.obsts(3-obst).par((5-2*obst)/4), 1, ...
            length(par.obsts(obst).colltau) ) ).^2, 1) ) ).*phi(par.obsts(obst).colltau, i);
%             exp(1i*par.k*sqrt( sum( (par.obsts(obst).par(par.obsts(obst).colltau) - repmat(par.obsts(3-obst).par((2*obst-1)/4), 1, ...
%             length(par.obsts(obst).colltau) ) ).^2, 1) ) ).*phi(par.obsts(obst).colltau, i);
%             length(par.obsts(obst).colltau) ) ).^2, 1) ) );
    end
end
sel = @(x,s) x(s);
phasFirRefl = zeros(length(par.obsts(1).colltau),4);
if nrPh > 2
    for obst = 1:length(par.obsts)
        cot = par.obsts(obst).colltau;
%             A(length(cot)*(obst-1)+1:length(cot)*obst, (length(par.obsts)*2+obst-1)*fN+i) = ...
        for coi = 1:length(cot)
            pt = par.obsts(obst).par(cot(coi)); % Find minimal distance from source to pt via other obst
%             [tauOth, phasi] = fminbnd( @(tau) sel(par.obsts(3-obst).par(tau),1) + sum( (pt - par.obsts(3-obst).par(tau) ).^2, 1), (2*obst-2)/4, (2*obst-0)/4);
            [tauOth, phasi] = fminbnd( @(tau) sel(par.obsts(3-obst).par(tau),1) + sum( (pt - par.obsts(3-obst).par(tau) ).^2, 1), (4-2*obst)/4, (6-2*obst)/4);
            phasFirRefl(coi,obst) = phasi;
            phasFirRefl(coi,obst+2) = tauOth;
%             phasi = par.obsts(3-obst).par(tauOth) + ; 
            % phasi is x-coord at other obst + dist from there to pt
            for i = 1:fN
                A(length(cot)*(obst-1)+coi, (length(par.obsts)*2+obst-1)*fN+i) = exp(1i*par.k*phasi).*phi(cot(coi), i);
            end
        end
    end
end
proj = A\c1;
figure; plot([real(proj) imag(proj)]); title('proj');

