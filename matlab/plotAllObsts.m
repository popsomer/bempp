% Plot all obstacles.
%% Initialising
clearvars
close all
clc
format longe
set(0,'DefaultFigureWindowStyle','docked');


%% All obstacles
mnr = 200;
ts = linspace(0,1,mnr)';
showTau = (0:5)/6;
labels = {'Circle', 'Ellipse', 'Near-inclusion', 'Almost convex', ...
    'Two circles', 'Three ellipses', 'Circle and near-inclusion', 'Nonconvex polygon'};
las = {'\xi=4%, N/k=6', '\xi=4%, N/k=6', '\xi=1%, N/k=10', '\xi=4%, N/k=7', ...
    '\xi=0.3%, N/k/2=10', '\xi=0.3%, N/k/3=10', '\xi=0.1%, N/k/2=10', '\xi=4%, N/k=9'};
hormove = [0.3, 0.05,0, -0.2  , -1, -0.2, -0.3, 0];
vertmove = [0, -0.1, 0.4, -0.1  , 0.5 , 1.1, 0.2, 0];

cols = 'rgbmkcy';
marks = {'v', '+', 'o', 'x', '*', 'h', 'd'};
figure; 
for oi = 1:8
    par = getObst(oi);
    subplot(2,4,oi);
    if isfield(par,'obsts')
        parametr = [];
        for moi = 1:length(par.obsts)
            topl = par.obsts(moi).par(ts');
            plot(topl(1,:), topl(2,:));
            hold on;
            h = text(mean(topl(1,:)) - 0.5*(max(topl(1,:))-min(topl(1,:)))+0.1, mean(topl(2,:)), ['Obst. ' num2str(moi)]);
            set(h,'FontSize',10);
            parametr = [parametr, topl];
            topl = par.obsts(moi).par(showTau);
            if (oi == 6) && (moi == 3)
                pt = par.obsts(moi).par(0.5);
                plot(pt(1), pt(2), 'ro','MarkerFaceColor','r', 'MarkerSize',5);
                h = text(pt(1)-0.3, pt(2)+0.1, '$\mathbf{p}$', 'interpreter', 'latex');
                set(h,'FontSize',10, 'color', 'r');
            end
        end
    else
        topl = par.par(showTau);
        parametr = par.par(ts');
        plot(parametr(1,:), parametr(2,:));
        if oi == 3
            hold on;
            pt = par.par(0.35);
            plot(pt(1), pt(2), 'ro','MarkerFaceColor','r', 'MarkerSize',5);
            h = text(pt(1)+0.05, pt(2)+0.04, '$\mathbf{q}$', 'interpreter', 'latex');
            set(h,'FontSize',10, 'color', 'r');
            
            pt = par.par(0.5);
            plot(pt(1), pt(2), 'ro','MarkerFaceColor','r', 'MarkerSize',5);
            h = text(pt(1)+0.02, pt(2)-0.02, '$\mathbf{r}$', 'interpreter', 'latex');
            set(h,'FontSize',10, 'color', 'r');
            
            pt = par.par(0.65);
            plot(pt(1), pt(2), 'ro','MarkerFaceColor','r', 'MarkerSize',5);
            h = text(pt(1), pt(2)-0.07, '$\mathbf{s}$', 'interpreter', 'latex');
            set(h,'FontSize',10, 'color', 'r');
        end
    end
    siz = max(parametr, [],2)-min(parametr,[],2);
    mid = (max(parametr, [],2) +min(parametr,[],2) )/2;
    hold on;
    rad = norm(siz);
    if abs(-1*besselh(0,1, 100*sqrt( (0.7 - 1).^2 + (pi -1).^2) )/par.bc(100, [0.7; pi]) -1 ) < 1e-12
        % Point source
        quiver(1 + 0.04*rad*cos(linspace(0,2*pi,6)), 1 + 0.04*rad*sin(linspace(0,2*pi,6)), ...
            0.1*rad*cos(linspace(0,2*pi,6)), 0.1*rad*sin(linspace(0,2*pi,6)), 'm', 'Autoscale','off');
    elseif abs(-1*exp(1i*222*([exp(1); -0.11]')*[cos(0); sin(0)])/par.bc(222, [exp(1); -0.11]) -1 ) < 1e-12
        % Incident plane wave
        quiver(mid(1) -0.8*siz(1)*ones(1,6), mid(2)+rad*linspace(-0.3,0.3,6), 0.1*rad*ones(1,6), zeros(1,6), 'm','Autoscale','off');
    else % superpos 3 plane waves
        quiver(mid(1)+0.8*rad*ones(1,6), mid(2)+siz(2)*linspace(-0.3,0.3,6), -0.3*rad*ones(1,6), zeros(1,6), 'm','Autoscale','off');
        quiver(mid(1)+ rad*0.7*cos(2*pi/3) +rad*linspace(-0.3,0.3,6)*cos(pi/6), mid(2) +rad*0.7*sin(2*pi/3) +rad*linspace(-0.3,0.3,6)*sin(pi/6), ...
            0.1*rad*cos(-pi/3)*ones(1,6), 0.1*rad*sin(-pi/3)*ones(1,6), 'm','Autoscale','off')
        quiver(mid(1) + rad*0.7*cos(-2*pi/3) +rad*linspace(-0.3,0.3,6)*cos(-pi/6), mid(2) +rad*0.7*sin(-2*pi/3) +rad*linspace(-0.3,0.3,6)*sin(-pi/6), ...
            0.2*rad*cos(pi/3)*ones(1,6), 0.2*rad*sin(pi/3)*ones(1,6), 'm','Autoscale','off')
    end
    title({labels{oi}; ['\rm{' las{oi} '}']});
    axis equal;
end
