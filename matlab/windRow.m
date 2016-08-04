% Sort, periodize and join the bounds and then construct a row of the collocation matrix using windows.
% Parts where the window == 1 can only become larger by joining windows. multElWind and multKernelWind are not mutually exclusive:
% if they are both zero, there is a block window.
% Input
%   i              - Row number
%   par            - The structure containing k,par and so on
%   bounds         - Each column gives the lambda, rho, l-lambda and rho-r of the window, but these could overlap or be too small.
%   [multElWind    - nonzero when multiplying the original matrix elements with the window function evaluated on the respective support, default 1]
%   [multKernelWind- nonzero when multiplying the kernel with the window function, default 0]
%   [collx         - Collocation point for when computing row in coupling matrix]
%   [normRois      - rois(i,:) divided by the threshold, in case windows are computed using correlation 'distances']
%   [ftz           - the factor below the threshold below which the window is zero, ftz=1 means block window, when correlation 'distances']
%   [roit          - the scaling between rois and A2, when correlation 'distances']
% Output
%   row            - a row of the collocation matrix
function row = windRow(i,par,bounds,multElWind, multKernelWind, collx,normRois, ftz,roit)
if ~exist('multElWind','var') || isempty(multElWind)
    multElWind = 1;
end
if ~exist('multKernelWind','var') || isempty(multKernelWind)
    multKernelWind = 0;
end
if ~exist('collx','var')
    collx = [];
end
windGreen = @(t) ones(size(t)); % The modification to the Green's function

%% Compute windows using correlation 'distances'
if exist('normRois','var') && ~isempty(normRois)
    getjr = @(t) arrayfun(@(tti) find(min(abs(roit-tti)) == abs(roit-tti),1,'first'), t);
    if multKernelWind
        windGreen = @(t) chi(abs(normRois(getjr(t))), ftz, 1,Inf,Inf,0);
    end
    if ~isempty(collx)
        row = collRowQBF('error', par, 1, par.N, windGreen, collx);
    else
        row = collRowQBF(i, par, 1, par.N, windGreen);
    end
    if multElWind
        row = chi(abs(normRois(getjr(par.colltau))), ftz, 1,Inf,Inf,0).*row;
    end
    return
end

%% Delete windows that are too small and then sort, periodize and join the bounds
saveBounds = bounds;

for nrWind = size(bounds,2):-1:1
    if bounds(2,nrWind) < bounds(1,nrWind)
        error('Did not expect for example [0.9; 0.1;..] = [0;0.1;..] & [0.9;1;..] at this point');
    elseif bounds(2,nrWind)-bounds(1,nrWind) < 3/par.N
        % The window is too small (maybe because the corresponding correlations are just above threshold): this would impede a Cinf window.
        if (isempty(collx) && (bounds(1,nrWind) < par.colltau(i) ) && (bounds(2,nrWind) > par.colltau(i) ) )
            warning('windRow:removeGreen', 'Removing the Green singularity because 3 points are not enough to capture it.');
        end
        bounds(:,nrWind) = [];
    end
end

if ~isempty(saveBounds)
    if isempty(bounds) % Take the largest window and extend it to include the closest colltau's.
        [siz, idxBiggest] = max(saveBounds(2,:)-saveBounds(1,:));
        if siz == 0
            error('Input bounds to windRow should not have end-T equal to start-T');
        end
        lambda = saveBounds(1,idxBiggest);
        rho = saveBounds(2,idxBiggest);
        cj1 = par.colltau(find(par.colltau <= lambda, 1, 'last'));
        if lambda < 0
            cj1 = par.colltau(find(par.colltau <= lambda+1, 1, 'last'))-1;
        elseif lambda > par.colltau(end)
            cj1 = par.colltau(find(par.colltau <= lambda-1, 1, 'last'))+1;
        end
        cj2 = par.colltau(find(par.colltau >= rho, 1, 'first'));
        if rho < 0
            cj2 = par.colltau(find(par.colltau >= rho+1, 1, 'first'))-1;
        elseif rho > par.colltau(end)
            cj2 = par.colltau(find(par.colltau >= rho-1, 1, 'first'))+1;
        end
        bounds = [cj1-1/par.N; cj2+1/par.N; (cj2-cj1+2/par.N)/siz*saveBounds(3:4,idxBiggest)];
    end
    bounds = [(bounds(1:2,:)-repmat(fix(bounds(1,:)),2,1) ); bounds(3:4,:)]; % Put the starts in [0,1]
    [~,I] = sort(bounds(1,:));
    bounds = bounds(:,I);
    
    if (bounds(1,1) <= 0) && (1+bounds(1,1) <= bounds(end,end) )
        % Take the negative variant
        bounds = [[bounds(1,end)-1; bounds(2:end,1)] bounds(:,2:end-1)];
    end
end

% 'bounds' should be sorted on the first row now. Check whether each next window can be joined with the previous one.
% We do this without a for-loop since the required indices may change when three windows must be joined
ol = 1;
while ol <= size(bounds,2)-1 % While there still is a next window to possibly join 'ol' with
    if(bounds(1,ol+1) - bounds(2,ol) <= 5/par.N +3*eps) % Windows overlapping or too close
        % b(1,ol+1) cannot be smaller than b(1,ol) and the periodic case is handled later.
        newbu = [bounds(1,ol); max(bounds(2,ol), bounds(2,ol+1))];
        % Make sure that newb(1,ol) +newb(3,ol) == min(b(1,ol)+b(3,ol), b(1,ol+1)+-b(3,ol+1)) :
        newbd = [min(bounds(1,ol)+bounds(3,ol), bounds(1,ol+1)+bounds(3,ol+1) ) - newbu(1); ...
            newbu(2)-max(bounds(2,ol)-bounds(4,ol), bounds(2,ol+1)-bounds(4,ol+1) )];
        % : makes sure that newb(2,ol) - newb(4,ol) == max(b(2,ol)-b(4,ol), b(2,ol+1)-b(4,ol+1)) so that the constant
        % part of the window can only be increased, though it should not be almost all of the window.
        bounds(:,ol) = [newbu; newbd];
        bounds(:,ol+1) = [];
        % no need to increase ol because size(bounds,2) decreased
    else
        ol = ol + 1;
    end
end

ol = 1;
while ol <= size(bounds,2)
    if(bounds(2,ol) - bounds(4,ol) <= bounds(1,ol)+bounds(3,ol) +9*eps) % Again remove windows which are too small
        % Windows can be too small when only one pt in roit above threshold or for artificial windows around tau = 0 or 1 due to I(ini)+1
        if isempty(collx) && (bounds(1,ol) < par.colltau(i)) && (bounds(2,ol) > par.colltau(i))
            warning('Removing window which contains Green singularity because it is too small to capture it.');
        end
        bounds(:,ol) = [];
    else
        ol = ol+1;
    end
end

% Also delete periodically, could be multiple times.
while (~isempty(bounds)) && (1+bounds(1,1) <= bounds(2,end) + 5/par.N + 3*eps)
    newbu = [min(bounds(1,end)-1,bounds(1,1)); max(bounds(2,1), bounds(2,end)-1)];
    newbd = [min(bounds(1,1)+bounds(3,1), bounds(1,end)-1+bounds(3,end) ) - newbu(1); ...
        newbu(2)-max(bounds(2,1)-bounds(4,1), bounds(2,end)-1-bounds(4,end) )];
    bounds(:,1) = [newbu; newbd];
    bounds(:,end) = [];
    % bounds(2,1) may now be larger than bounds(1,2) so again loop to join
    ol = 1;
    while ol <= size(bounds,2)-1
        if(bounds(1,ol+1) - bounds(2,ol) <= 5/par.N +3*eps)
            newbu = [bounds(1,ol); max(bounds(2,ol), bounds(2,ol+1))];
            newbd = [min(bounds(1,ol)+bounds(3,ol), bounds(1,ol+1)+bounds(3,ol+1) ) - newbu(1); ...
                newbu(2)-max(bounds(2,ol)-bounds(4,ol), bounds(2,ol+1)-bounds(4,ol+1) )];
            bounds(:,ol) = [newbu; newbd];
            bounds(:,ol+1) = [];
        else
            ol = ol + 1;
        end
    end
end

%% Construct the row
% Also take colltau into account: when [-1e-12; 0.9999;...] then all collocation points fall into the window, so we make the full row.
if isempty(bounds) || (max(bounds(2,:)-bounds(1,:)) > par.colltau(end))
    if exist('collx','var') && ~isempty(collx)
        row = collRowQBF('error',par,1,par.N,[],collx);
    else
        row = collRowQBF(i,par);
    end
    return
end
row = zeros(1,par.N);
for nrWind = 1:size(bounds,2)
    lambda = bounds(1,nrWind);
    rho = bounds(2,nrWind);
    l = lambda + bounds(3,nrWind);
    r = rho - bounds(4,nrWind);
    if lambda > rho
        error('Expected rho >= lambda');
    elseif l > r
        error('Constant 1 part of window wrong in bounds');
    elseif (nrWind >= 2) && (lambda < bounds(2,nrWind-1))
        error('Bounds were not correctly joined');
    end
    
    scp = @(t) chi(t,lambda,l, r,rho,1); % Smooth Cutoff Periodic
    % We could test whether scp(par.colltau) has min 0 and max 1 but this might not be the case for a small window or l close to r.
    if multKernelWind
        windGreen = scp;
    end
    [j1, j2, noSplit] = bounds2Ind(par.colltau,lambda,rho);
    if noSplit % j1 = j2 or j2-1 should have have been removed when windows are too small, so now colltau(j1) < lambda < rho < colltau(j2)
        if  exist('collx','var') && ~isempty(collx)
            row(j1:j2) = collRowQBF('error', par, j1, j2, windGreen, collx);
        else
            row(j1:j2) = collRowQBF(i, par, j1, j2, windGreen);
        end
        if multElWind
            row(j1:j2) = scp(par.colltau(j1:j2)).*row(j1:j2);
        end
    else % A window like (-0.1,0.2) has to be split.
        if  exist('collx','var') && ~isempty(collx)
            row(j1:par.N) = collRowQBF('error',par,j1,par.N,windGreen,collx);
            row(1:j2) = collRowQBF('error',par,1,j2,windGreen,collx);
        else
            row(j1:par.N) = collRowQBF(i,par,j1,par.N,windGreen);
            row(1:j2) = collRowQBF(i,par,1,j2,windGreen);
        end
        if multElWind
            row(j1:par.N) = scp(par.colltau(j1:par.N)).*row(j1:par.N);
            row(1:j2) = scp(par.colltau(1:j2)).*row(1:j2);
        end
    end
end

end
