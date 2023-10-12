% Computes n-dimensional vector fields for potential system with n
% reference states, starting in the fixed frame

% PotentialField2(2,'steps',20,'startPoint',[1,1],'npoints',1000,'step_size',1e-3,'save',true)

function PotentialField2(varargin)

p = inputParser;
p.CaseSensitive = 1;

addRequired(p,'n',@(X) X >= 2)
addParameter(p,'c',sqrt(2),@(X) X > 0); % distance between each pair of states
addParameter(p,'steps',10, @(X) X >= 1); % number of divisions of c
addParameter(p,'extraStates',[]); % matrix of extra reference states
addParameter(p,'startPoint',[]); % number of divisions of c
addParameter(p,'npoints',1000); % number of integration points
addParameter(p,'step_size',1e-3); % integration step size
addParameter(p,'save',false); % whether to save image or not

parse(p,varargin{:})

n = p.Results.n; % dimension of space, initial number of reference states
c = p.Results.c; % distance between reference states
steps = p.Results.steps; % number of divisions of c
npoints = p.Results.npoints; % number of integration points
step_size = p.Results.step_size; % number of integration points
save = p.Results.save; % whether to save image or not

dc = c / steps; 

States = eye(n) * c / sqrt(2); % coordinates of reference states
if ~isempty(p.Results.extraStates) && size(p.Results.extraStates,2) == n
    States = cat(1,States,p.Results.extraStates);
    m = size(States,1); % new number of reference states
else
    m = n;
end

char = strlength(num2str(m)); % number of characters in the number m of states
str = '[x1'; % create m-dimensional mesh
for i = 2:n
    fmt = sprintf('%%s,x%%0.%sd',num2str(char));
    str = sprintf(fmt,str,i);
end
str = sprintf('%s] = ndgrid(0:dc:c);',str);
eval(str);

x = cell(m,1);
[x{:}] = deal(zeros(size(x1)));
for i = 1:n % transfer X's into a cell array
    eval(sprintf('x{%d} = x%d;',i,i))
end

e = cell(m,n); % basis vectors in relative frame
X = cell(m,1); % relative frame coordinates
for i = 1:m
    for j = 1:n
        e{i,j} = x{j} - States(i,j); % jth component of ith basis vector 
    end
    X{i} = sqrt(sum(cat(n + 1,e{i,:}) .^ 2,n + 1));
    for j = 1:n
        e{i,j} = e{i,j} ./ X{i}; % normalize basis vectors
    end
end

V = cell(m,1); % relative frame velocities
[V{:}] = deal(zeros(size(x1)));
for i = 1:m % define ith velocity component
    for j = 1:m
        if i < j
            V{i} = V{i} + (-1) ^ (i + j - 1) * X{j};
        elseif i > j
            V{i} = V{i} + (-1) ^ (i + j) * X{j};
        end
    end
end

v = cell(n,1); % fixed frame velocities
[v{:}] = deal(zeros(size(x1)));
for i = 1:n
    for j = 1:m
        v{i} = v{i} + V{j} .* e{j,i}; % ith fixed frame velocity due to jth relative frame velocity
    end
end

Pairs = nchoosek(1:m,2); % Pairs of reference states

H = sum(cat(n + 1,X{:}) .^ 2,n + 1) / 2; %energy function

Phi = arrayfun(@(I) X{Pairs(I,1)} .* X{Pairs(I,2)},1:nchoosek(m,2),'UniformOutput',false);
Phi = sum(cat(n + 1,Phi{:}),n + 1); % potential function

Psi = arrayfun(@(I) (X{Pairs(I,1)} .^ 2 - X{Pairs(I,2)} .^ 2) / 2,1:nchoosek(m,2),'UniformOutput',false);
Psi = sum(cat(n + 1,Psi{:}),n + 1); % conjugate potential


% integral curve
if isempty(p.Results.startPoint) || size(p.Results.startPoint,2) ~= n
    c0 = States(1,:) + abs(norminv(rand(1,n),0,0.001)); % initial point with small positive normal error
    fprintf('Starting point is: [')
    fprintf('%s, ',c0(1:end - 1))
    fprintf('%s]\n,',c0(end))
else
    c0 = p.Results.startPoint; % supplied initial point
    fprintf('Starting point is: [')
    fprintf('%f, ',c0(1:end - 1))
    fprintf('%f]\n,',c0(end))
end

E_new = arrayfun(@(I) (c0 - States(I,:)) / sqrt(sum((c0 - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
X_new = arrayfun(@(I) sqrt(sum((c0 - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
phi_new = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(m,2),'UniformOutput',false);
phi_new = sum(cat(1,phi_new{:}),1); % potential function
h_new = arrayfun(@(I) X_new{I} .^ 2 / 2,1:m,'UniformOutput',false);
h_new = sum(cat(1,h_new{:}),1); % energy function
V_new = cell(m,1);
[V_new{:}] = deal(zeros(1,n));
for i = 1:m % define ith velocity component
    for j = 1:m
        if i < j
            V_new{i} = V_new{i} + (-1) ^ (i + j - 1) * X_new{j} * E_new{i}; % velocity along ith basis vector due to jth distance
        elseif i > j
            V_new{i} = V_new{i} + (-1) ^ (i + j) * X_new{j} * E_new{i};
        end
    end
end

gamma = zeros(n,npoints); % initialize matrix of sample points
gamma(:,1) = c0'; % first point
phi = zeros(1,npoints);
phi(1) = phi_new;
h = zeros(1,npoints);
h(1) = h_new;
for k = 2:npoints
    c_new = gamma(:,k - 1)'; % last point
    V_old = V_new;
    for i = 1:m
        c_new = c_new + V_old{i} * step_size; % update position
    end
    gamma(:,k) = c_new';
    E_new = arrayfun(@(I) (c_new - States(I,:)) / sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
    X_new = arrayfun(@(I) sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
    phi_new = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(m,2),'UniformOutput',false);
    phi_new = sum(cat(1,phi_new{:}),1); % potential function
    phi(k) = phi_new;
    h_new = arrayfun(@(I) X_new{I} .^ 2 / 2,1:m,'UniformOutput',false);
    h_new = sum(cat(1,h_new{:}),1); % energy function
    h(k) = h_new;
    V_new = cell(m,1);
    [V_new{:}] = deal(zeros(1,n));
    for i = 1:m % define ith velocity component
        for j = 1:m
            if i < j
                V_new{i} = V_new{i} + (-1) ^ (i + j - 1) * X_new{j} * E_new{i}; % velocity along ith basis vector due to jth distance
            elseif i > j
                V_new{i} = V_new{i} + (-1) ^ (i + j) * X_new{j} * E_new{i};
            end
        end
    end
end

% reverse trajectory
E_new = arrayfun(@(I) (c0 - States(I,:)) / sqrt(sum((c0 - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
X_new = arrayfun(@(I) sqrt(sum((c0 - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
phi_new_rev = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(m,2),'UniformOutput',false);
phi_new_rev = sum(cat(1,phi_new_rev{:}),1); % potential function
h_new_rev = arrayfun(@(I) X_new{I} .^ 2 / 2,1:m,'UniformOutput',false);
h_new_rev = sum(cat(1,h_new_rev{:}),1); % energy function
V_new = cell(m,1);
[V_new{:}] = deal(zeros(1,n));
for i = 1:m % define ith velocity component
    for j = 1:m
        if i < j
            V_new{i} = V_new{i} + (-1) ^ (i + j - 1) * X_new{j} * E_new{i}; % velocity along ith basis vector due to jth distance
        elseif i > j
            V_new{i} = V_new{i} + (-1) ^ (i + j) * X_new{j} * E_new{i};
        end
    end
end

gamma_rev = zeros(n,npoints); % initialize matrix of sample points
gamma_rev(:,1) = c0'; % first point
phi_rev = zeros(1,npoints);
phi_rev(1) = phi_new_rev;
h_rev = zeros(1,npoints);
h_rev(1) = h_new_rev;
for k = 2:npoints
    c_new = gamma_rev(:,k - 1)'; % last point
    V_old = V_new;
    for i = 1:m
        c_new = c_new - V_old{i} * step_size; % update position
    end
    gamma_rev(:,k) = c_new';
    E_new = arrayfun(@(I) (c_new - States(I,:)) / sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
    X_new = arrayfun(@(I) sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
    phi_new_rev = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(m,2),'UniformOutput',false);
    phi_new_rev = sum(cat(1,phi_new_rev{:}),1); % potential function
    phi_rev(k) = phi_new_rev;
    h_new_rev = arrayfun(@(I) X_new{I} .^ 2 / 2,1:m,'UniformOutput',false);
    h_new_rev = sum(cat(1,h_new_rev{:}),1); % energy function
    h_rev(k) = h_new_rev;
    V_new = cell(m,1);
    [V_new{:}] = deal(zeros(1,n));
    for i = 1:m % define ith velocity component
        for j = 1:m
            if i < j
                V_new{i} = V_new{i} + (-1) ^ (i + j - 1) * X_new{j} * E_new{i}; % velocity along ith basis vector due to jth distance
            elseif i > j
                V_new{i} = V_new{i} + (-1) ^ (i + j) * X_new{j} * E_new{i};
            end
        end
    end
end

gamma = [fliplr(gamma_rev(:,2:end)),gamma];
phi = [fliplr(phi_rev(:,2:end)),phi];
h = [fliplr(h_rev(:,2:end)),h];

if n == 2
    if ishandle(1)
        set(0,'CurrentFigure',1)
        cla
    else
        figure(1)
    end
    
    surf(x{1},x{2},H,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.3);
    view(15,45)
    hold on
    % plot3(gamma(1,:),gamma(2,:),c ^ 2 / 2 * ones(size(gamma,2)),'LineWidth',2,'Color','k')
    plot3(gamma(1,:),gamma(2,:),h,'LineWidth',2,'Color','k')
    quiver3(x{1},x{2},H,v{1},v{2},zeros(size(H)),'Color','k','LineWidth',1);
    axis([0,c,0,c,0,Inf])
    % daspect([1,1,1])
    
    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.ZLabel.String = '$$H$$';
    ax.ZLabel.Rotation = 0;
    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
    set(get(gca,'XLabel'),'Interpreter','latex')
    set(get(gca,'YLabel'),'Interpreter','latex')
    set(get(gca,'ZLabel'),'Interpreter','latex')
    hold off

    if save
        exportgraphics(gcf,sprintf('%s/Trajectory n = %d, H surface, start = (%0.6f,%0.6f).png',pwd,n,c0(1),c0(2)),'Resolution',300)
    end

    if ishandle(2)
        set(0,'CurrentFigure',2)
        cla
    else
        figure(2)
    end
    
    surf(x{1},x{2},H,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.3);
    view(0,90) % top-down view of energy surface
    cb = colorbar;
    set(cb,'Limits',[0,max(H(:))])
    cb.Label.String = '$$H$$';
    cb.Label.Rotation = 0;
    cb.Label.Interpreter = 'latex';
    set(cb,'FontSize',16)
    hold on
    % plot3(gamma(1,:),gamma(2,:),c ^ 2 / 2 * ones(size(gamma,2)),'LineWidth',2,'Color','k')
    plot3(gamma(1,:),gamma(2,:),h,'LineWidth',2,'Color','k')
    scatter3(States(:,1),States(:,2),zeros(size(States,1),1),200,...
        'filled','r','o','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
    quiver3(x{1},x{2},H,v{1},v{2},zeros(size(H)),'Color','k','LineWidth',1);
    axis([0,c,0,c,0,Inf])
    % daspect([1,1,1])
    
    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.YLabel.Rotation = 0;
    ax.ZLabel.String = '$$H$$';
    ax.ZLabel.Rotation = 0;
    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
    set(get(gca,'XLabel'),'Interpreter','latex')
    set(get(gca,'YLabel'),'Interpreter','latex')
    set(get(gca,'ZLabel'),'Interpreter','latex')
    hold off

    if save
        exportgraphics(gcf,sprintf('%s/Trajectory n = %d, H overlaid, start = (%0.6f,%0.6f).png',pwd,n,c0(1),c0(2)),'Resolution',300)
    end

    if ishandle(3)
        set(0,'CurrentFigure',3)
        clf
    else
        figure(3)
    end
    % plot trajectory and vector field on potential energy surface
    surf(x{1},x{2},Phi,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.3);
    view(15,45)
    hold on
    plot3(gamma(1,:),gamma(2,:),phi,'LineWidth',2,'Color','k')
    quiver3(x{1},x{2},Phi,v{1},v{2},-2 * Psi,'Color','k','LineWidth',1);
    axis([0,c,0,c,0,Inf])
    % daspect([1,1,1])
    
    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.ZLabel.String = '$$\Phi$$';
    ax.ZLabel.Rotation = 0;
    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
    set(get(gca,'XLabel'),'Interpreter','latex')
    set(get(gca,'YLabel'),'Interpreter','latex')
    set(get(gca,'ZLabel'),'Interpreter','latex')
    hold off

    if save
        exportgraphics(gcf,sprintf('%s/Trajectory n = %d, Phi surface, start = (%0.6f,%0.6f).png',pwd,n,c0(1),c0(2)),'Resolution',300)
    end

    if ishandle(4)
        set(0,'CurrentFigure',4)
        clf
    else
        figure(4)
    end
    % plot trajectory and vector field on potential energy surface
    surf(x{1},x{2},Phi,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.3);
    cb = colorbar;
    set(cb,'Limits',[0,max(Phi(:))])
    cb.Label.String = '$$\Phi$$';
    cb.Label.Rotation = 0;
    cb.Label.Interpreter = 'latex';
    set(cb,'FontSize',16)
    view(0,90) % top-down view of potential surface
    hold on
    plot3(gamma(1,:),gamma(2,:),phi,'LineWidth',2,'Color','k')
    scatter3(States(:,1),States(:,2),zeros(size(States,1),1),200,...
        'filled','r','o','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
    quiver3(x{1},x{2},Phi,v{1},v{2},-2 * Psi,'Color','k','LineWidth',1);
    axis([0,c,0,c,0,Inf])
    % daspect([1,1,1])
    
    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.YLabel.Rotation = 0;
    ax.ZLabel.String = '$$\Phi$$';
    ax.ZLabel.Rotation = 0;
    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
    set(get(gca,'XLabel'),'Interpreter','latex')
    set(get(gca,'YLabel'),'Interpreter','latex')
    set(get(gca,'ZLabel'),'Interpreter','latex')
    hold off

    if save
        exportgraphics(gcf,sprintf('%s/Trajectory n = %d, Phi overlaid, start = (%0.6f,%0.6f).png',pwd,n,c0(1),c0(2)),'Resolution',300)
    end

end

