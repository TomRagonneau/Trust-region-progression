% Generate a sequence of figures that represents the course of a generic
% trust-region optimization method.
%
% The execution of this script requires MATLAB R2020a or later.
%
% Copyright 2020 Tom M. Ragonneau.
%
close all;
set(gcf,'Position',get(0,'Screensize'),'Visible','off');
set(0,'defaultTextInterpreter','latex');

% Define some constants related to the trust-region algorithm.
radiusmax = 3;  % maximal trust-region radius
radiusmin = .05;  % minimal trust-region radius (stop the execution)
radius = .75;  % initial trust-region radius
xopt = 0;  % first coordinate of the initial guess
yopt = -2;  % second coordinate of the initial guess
xmin = -4.5;  % minimal value of the first coordinate
xmax = 4.5;  % maximal value of the first coordinate
ymin = -3;  % minimal value of the second coordinate
ymax = 3;  % maximal value of the second coordinate
kmax = 15;  % maximal number of outer-loop iteration

dpi = 200;  % desired resolution of the images (in dots per inch)
fontsize = 15;  % font size of the legends, axis labels and colorbar labels
step = .01;  % step of the mesh of computation
colors = get(gca,'colororder');  % default MATLAB colors

% The figures are saved in the folder `saving`, and will be saved as
% <saving>/<fname><no>.png.
saving = 'imgs';
fname = 'tr';
if ~exist(saving, 'dir')
    mkdir(saving)
end

% Generate the mesh of computation and evaluate the objective function.
[X,Y] = meshgrid(xmin:step:xmax,ymin:step:ymax);
xyopt = intersect(find(abs(X-xopt) < step/2),find(abs(Y-yopt) < step/2));
if length(xyopt) ~= 1
    % There should exists positive integers k1 and k2 such that
    % xopt = xmin+k1*step and yopt=ymin+k2*step, so that the initial guess
    % belongs to the mesh.
    ME = MException('TrustRegion:InitializationError','Ensure that xopt and yopt belong to X and Y.');
    throw(ME);
end
Z = objective(X,Y);

% Computes the index of the minimum of the objective function on the mesh.
% Note that its coordinates can be obtained by [X(xysol); Y(xysol)]. The
% index correspondst to the absolute index of the matrices, columnwise.
[~,xysol] = min(Z,[],'all','linear');
Zmin = min(Z,[],'all');
Zmax = max(Z,[],'all');

% Generate the initial figure, that represents the objective function and
% its minimum. The legend contains only the information related to the
% isopleth graph.
contour(X,Y,Z,40,'Linewidth',2);
colorbar('FontSize',fontsize,'TickLabelInterpreter','latex');
colormap pink;
caxis([Zmin,Zmax]);
axis equal;
hold on;
scatter(X(xysol),Y(xysol),250,'k','h','MarkerFaceColor',colors(7,:),'LineWidth',1.5);
xlim([xmin,xmax]);
ylim([ymin,ymax]);
hold off;
legend({'Objective function'},'Location','northeast','FontSize',fontsize,'Interpreter','latex');
set(gca,'Xticklabel',get(gca,'Xticklabel'),'FontSize',fontsize,'TickLabelInterpreter','latex');
exportgraphics(gcf,sprintf('%s/%s1.png',saving,fname),'Resolution',dpi);

% Start the trust-region iterations. The values of the model are computed
% only if the best point so far has been calculated at the previous
% iteration. The variables Zm and Zt will be initialized with the values of
% the model during the first iteration.
xyplus = xyopt;
Zm = [];
Zt = [];
for k = 1:kmax
    % Compute the Taylor expansion of the objective function around the
    % best point so far. The value of the objective function is set to NaN
    % in the trust region not to be drawn.
    if xyplus == xyopt
        Zm = model(X,Y,xopt,yopt);
        Zt = Z(:,:);
        Zmb = Zm(:,:);
        Zt((X-xopt).^2+(Y-yopt).^2 <= radius^2) = NaN;
        Zmb((X-xopt).^2+(Y-yopt).^2 > radius^2) = NaN;

        % Compute the index of the minimum reached by the model. Note that
        % the variable Zmb contains NaN at each point outside of the trust
        % region, so that the minimum is computed in the trust-region ball.
        [~,xyplus] = min(Zmb,[],'all','linear');
    end
    
    % Generate the figure that represents the current trust-region
    % subproblem, without indicating the minimum reached by the model.
    lgd = {'Objective function','Quadratic model','Solution~$\underline{x}^{\ast}$','Iterate~$\underline{x}_k$'};
    contour(X,Y,Zt,40,'Linewidth',2);
    colorbar('FontSize',fontsize,'TickLabelInterpreter','latex');
    colormap pink;
    caxis([Zmin,Zmax]);
    axis equal;
    hold on;
    contour(X,Y,Zm,40,'--','Linewidth',2);
    scatter(X(xysol),Y(xysol),250,'k','h','MarkerFaceColor',colors(7,:),'LineWidth',1.5);
    scatter(xopt,yopt,150,'k','o','MarkerFaceColor',colors(5,:),'LineWidth',1.5);
    rectangle('Position',[xopt-radius,yopt-radius,2*radius,2*radius],'Curvature',[1,1],'LineWidth',3);
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    legend(lgd,'Location','northeast','FontSize',fontsize,'Interpreter','latex');
    set(gca,'Xticklabel',get(gca,'Xticklabel'),'FontSize',fontsize,'TickLabelInterpreter','latex');
    exportgraphics(gcf,sprintf('%s/%s%.0f.png',saving,fname,2*k),'Resolution',dpi);
    
    % Generate the figure that represents the current trust-region
    % subproblem and indicate the minimum reached by the model. The legend
    % contains only the information related to the both isopleth graphs.
    % The legend is updated to insert the new point.
    lgd{end+1} = 'Trial point~$\underline{x}_k + \underline{d}_k$ with~$\|\underline{d}_k\|_2 \leq \Delta_k$';
    scatter(X(xyplus),Y(xyplus),250,'k','h','MarkerFaceColor',colors(2,:),'LineWidth',1.5);
    hold off;
    legend(lgd,'Location','northeast','FontSize',fontsize,'Interpreter','latex');
    set(gca,'Xticklabel',get(gca,'Xticklabel'),'FontSize',fontsize,'TickLabelInterpreter','latex');
    exportgraphics(gcf,sprintf('%s/%s%.0f.png',saving,fname,2*k+1),'Resolution',dpi);
    
    % Compute the trust-region updates. If the computed ratio is NaN, it
    % means that two consecutive values of the model were (almost) equal,
    % the ratio is therefore set to zero to reduce the trust-region radius.
    if radius <= radiusmin
        % The final trust-region radius has been reached, the progression
        % should stop.
        break;
    end
    ratio = (Z(xyopt)-Z(xyplus))/(Zm(xyopt)-Zm(xyplus));
    if isnan(ratio)
        ratio = 0;
    end
    fprintf('Iteration no %-2.0f: RADIUS = %.4f, RATIO = %.4f\n',k,radius,ratio)
    if ratio < 1/4
        radius = max(radiusmin,radius/3);
    elseif (ratio > 3/4) && ((X(xyopt)-X(xyplus))^2+(Y(xyopt)-Y(xyplus))^2 >= radius^2-2*step)
        % This condition is reached when a substancial true decrease has
        % been observed in the objective function and the computed trial
        % point is on the trust-region boundary.
        radius = min(radiusmax,2*radius);
    end
    if ratio >= 1e-1
        % Update the current best point so far.
        xyopt = xyplus;
        xopt = X(xyopt);
        yopt = Y(xyopt);
    end
end


function fx = objective(x,y)
%OBJECTIVE A two-dimensional objective function.
%   The considered function is f(x) = 10*(y^2-x^2)+12*sin(x*y)-2*x+x^4.
%
fx = 10*(y.^2-x.^2)+12*sin(x.*y)-2*x+x.^4;

return;
end


function mx = model(x,y,x0,y0)
%MODEL The model of the objective function around [x0; y0].
%   This function returns the second-order Taylor expansion of the
%   objective function.
%

% Zero-order term.
mx = objective(x0,y0);

% First-order term.
gx = 4*x0^3-20*x0+12*y0*cos(x0*y0)-2;
gy = 20*y0+12*x0*cos(x0*y0);
mx = mx+gx*(x-x0)+gy*(y-y0);

% Second-order term.
Hxx = 12*x0^2-12*y0^2*sin(x0*y0)-20;
Hxy = 12*(cos(x0*y0)-x0*y0*sin(x0*y0));
Hyy = 20-12*x0^2*sin(x0*y0);
mx = mx+(Hxx*(x-x0).^2+Hyy*(y-y0).^2+2*Hxy*(x-x0).*(y-y0))/2;

return;
end

