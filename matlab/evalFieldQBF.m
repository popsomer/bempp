% Evaluate the scattered field.
% Input
%   par		  - The structure containing k,N,par,gradnorm,corners,t,colltau
%   x		  - Point where to evaluate
%   sol		  - Solution vector
%   comprSol 	  - 1 if the solution is compressed (for possibly adding the phase), else 0
% Output
%   z		  - The scattered field
function z = evalFieldQBF(par,x,sol,comprSol)

if isfield(par,'obsts')
    par.dbf = par.obsts(1).dbf; % Assume the same degree of basis functions to avoid refactoring.
end

if par.dbf == 1
	qbf_w = [   0.01666666666667   0.26666666666667   0.43333333333333   0.26666666666667   0.01666666666667];
	step = 1/2; istep = 1/step;
    dbf = 1;
elseif par.dbf == 3  % Weights computed in iepack/basis/make_bfinfo_periodic_spline_equi.m using quad_sf with lo_d = sqrt(2)/16*[1 4 6 4 1]
	qbf_w = [  0.000022045855379   0.009876543209876   0.084744268077601  0.238447971781305   0.333818342151675  ...
		0.238447971781305  0.084744268077601   0.009876543209877   0.000022045855379];
	step = 4/(length(qbf_w)-1);
	istep = round(1/step);
    dbf = 3;
end
if isfield(par,'obsts') % Multiple scattering
	z = 0;
	for obst = 1:length(par.obsts)
		tau = ((0:par.obsts(obst).N*istep-1)*step -dbf)/par.obsts(obst).N;
		% Add points to the right using periodicity.
		tau2 = [tau tau(1:dbf*istep + 1)];
        
		kernelvals = 1i/4*besselh(0,1,par.k*sqrt(sum((repmat(x,1,size(tau2,2))-par.obsts(obst).par(tau2) ).^2, 1)));
		gnvals = par.obsts(obst).gradnorm(tau2);
		zt = 0;
		for i=1:par.obsts(obst).N
			zt = zt + sol(i+par.r(1,obst)-1) * sum(qbf_w .* kernelvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)) ...
				.* gnvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)));
		end
		z = z + zt/par.obsts(obst).N;
	end
elseif comprSol && isfield(par,'difase')
    fN = length(par.fco);
	tau = ((0:fN*istep-1)*step -dbf)/fN;
	% Add points to the right using periodicity.
	tau2 = [tau tau(1:dbf*istep + 1)];
	
	kernelvals = 1i/4*besselh(0,1,par.k*sqrt(sum((repmat(x,1,size(tau2,2))-par.par(tau2) ).^2, 1)));
	gnvals = par.gradnorm(tau2).*exp(2i*pi*par.k*interp1(par.Tase,par.difase,tau2, 'spline', 'extrap'));
	z = 0;
	for i=1:fN
		z = z + sol(i) * sum(qbf_w .* kernelvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)) ...
			.* gnvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)));
	end
	z = z/fN;
elseif comprSol && isfield(par,'phase')
    fN = length(par.fco);
	tau = ((0:fN*istep-1)*step -dbf)/fN;
	% Add points to the right using periodicity.
	tau2 = [tau tau(1:dbf*istep + 1)];
	
	kernelvals = 1i/4*besselh(0,1,par.k*sqrt(sum((repmat(x,1,size(tau2,2))-par.par(tau2) ).^2, 1)));
	gnvals = par.gradnorm(tau2).*exp(1i*par.k*par.phase(tau2));
	z = 0;
	for i=1:fN
		z = z + sol(i) * sum(qbf_w .* kernelvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)) ...
			.* gnvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)));
	end
	z = z/fN;
else
	tau = ((0:par.N*istep-1)*step -dbf)/par.N;
	% Add points to the right using periodicity.
	tau2 = [tau tau(1:dbf*istep + 1)];
	
	kernelvals = 1i/4*besselh(0,1,par.k*sqrt(sum((repmat(x,1,size(tau2,2))-par.par(tau2) ).^2, 1)));
	gnvals = par.gradnorm(tau2);
	z = 0;
	for i=1:par.N
		z = z + sol(i) * sum(qbf_w .* kernelvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)) ...
			.* gnvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)));
	end
	z = z/par.N;
end
