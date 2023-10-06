% Computes the difference between the n+1- and n-dimensional fields, i.e.,
% how an n-dimensional trajectory morphs into an n+1-dimensional trajectory

% PotentialField3(2,'steps',20,'startPoint',[1,1],'extraStates',[1,1],'npoints',1000,'step_size',1e-3,'magnification',100,'save',true)
% PotentialField3(2,'steps',20,'startPoint',[1,1],'extraStates',[1 + sin(15 * pi / 180),1 + sin(15 * pi / 180)],'npoints',1000,'step_size',1e-3,'magnification',100,'save','true')
% PotentialField3(2,'steps',20,'startPoint',[0.5,0.5],'extraStates',[1,1],'npoints',1000,'step_size',1e-3,'magnification',100,'save',true)

function PotentialField3(varargin)

p = inputParser;
p.CaseSensitive = 1;

addRequired(p,'n',@(X) X >= 2)
addParameter(p,'c',sqrt(2),@(X) X > 0); % distance between each pair of states
addParameter(p,'steps',10, @(X) X >= 1); % number of divisions of c
addParameter(p,'extraStates',[]); % matrix of extra reference states
addParameter(p,'startPoint',[]); % number of divisions of c
addParameter(p,'npoints',1000); % number of integration points
addParameter(p,'step_size',1e-3); % integration step size
addParameter(p,'magnification',10); % magnification factor
addParameter(p,'save',false); % whether to save image or not

parse(p,varargin{:})

n = p.Results.n; % dimension of space, initial number of reference states
c = p.Results.c; % distance between reference states
steps = p.Results.steps; % number of divisions of c
npoints = p.Results.npoints; % number of integration points
step_size = p.Results.step_size; % number of integration points
magnification = p.Results.magnification; % magnification factor
save = p.Results.save;

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

Y = cell(n,1); % relative frame coordinates
for i = 1:n % only for n-state system
    for j = 1:n
        e{i,j} = x{j} - States(i,j); % jth component of ith basis vector 
    end
    Y{i} = sqrt(sum(cat(n + 1,e{i,:}) .^ 2,n + 1));
    for j = 1:n
        e{i,j} = e{i,j} ./ Y{i}; % normalize basis vectors
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

W = cell(n,1); % relative frame velocities, in n-dimensions only
[W{:}] = deal(zeros(size(x1)));
for i = 1:n % define ith velocity component
    for j = 1:n
        if i < j
            W{i} = W{i} + (-1) ^ (i + j - 1) * X{j};
        elseif i > j
            W{i} = W{i} + (-1) ^ (i + j) * X{j};
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

w = cell(n,1); % fixed frame velocities, in only n dimensions
[w{:}] = deal(zeros(size(x1)));
for i = 1:n
    for j = 1:n
        w{i} = w{i} + W{j} .* e{j,i}; % ith fixed frame velocity due to jth relative frame velocity
    end
end

Pairs = nchoosek(1:m,2); % Pairs of reference states

% H = sum(cat(n + 1,X{:}) .^ 2,n + 1) / 2; %energy function
H = sum(cat(n + 1,Y{:}) .^ 2,n + 1) / 2; %energy function, n-dimensional

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

W_new = cell(n,1);
[W_new{:}] = deal(zeros(1,n));
for i = 1:n % define ith velocity component, in n dimensions only
    for j = 1:n
        if i < j
            W_new{i} = W_new{i} + (-1) ^ (i + j - 1) * X_new{j} * E_new{i}; % velocity along ith basis vector due to jth distance
        elseif i > j
            W_new{i} = W_new{i} + (-1) ^ (i + j) * X_new{j} * E_new{i};
        end
    end
end

gamma = zeros(n,npoints); % initialize matrix of sample points
gamma(:,1) = c0'; % first point

gamma_prime = zeros(n,npoints); % initialize matrix of sample points
gamma_prime(:,1) = c0'; % first point

for k = 2:npoints
    c_new = gamma(:,k - 1)'; % last point
    d_new = gamma(:,k - 1)'; % last point
    V_old = V_new;
    W_old = W_new;
    for i = 1:n
        c_new = c_new + W_old{i} * step_size; % update position
    end
    gamma(:,k) = c_new';
    for i = 1:m
        d_new = d_new + V_old{i} * step_size * magnification; % updated position following m-dimensional field from current point
    end
    gamma_prime(:,k) = d_new';
    E_new = arrayfun(@(I) (c_new - States(I,:)) / sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
    X_new = arrayfun(@(I) sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates

    W_new = cell(n,1);
    [W_new{:}] = deal(zeros(1,n));
    for i = 1:n % define ith velocity component
        for j = 1:n
            if i < j
                W_new{i} = W_new{i} + (-1) ^ (i + j - 1) * X_new{j} * E_new{i}; % velocity along ith basis vector due to jth distance
            elseif i > j
                W_new{i} = W_new{i} + (-1) ^ (i + j) * X_new{j} * E_new{i};
            end
        end
    end

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

% get reverse trajectory
E_new = arrayfun(@(I) (c0 - States(I,:)) / sqrt(sum((c0 - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
X_new = arrayfun(@(I) sqrt(sum((c0 - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates

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

W_new = cell(n,1);
[W_new{:}] = deal(zeros(1,n));
for i = 1:n % define ith velocity component, in n dimensions only
    for j = 1:n
        if i < j
            W_new{i} = W_new{i} + (-1) ^ (i + j - 1) * X_new{j} * E_new{i}; % velocity along ith basis vector due to jth distance
        elseif i > j
            W_new{i} = W_new{i} + (-1) ^ (i + j) * X_new{j} * E_new{i};
        end
    end
end

gamma_rev = zeros(n,npoints); % initialize matrix of sample points
gamma_rev(:,1) = c0'; % first point

gamma_prime_rev = zeros(n,npoints); % initialize matrix of sample points
gamma_prime_rev(:,1) = c0'; % first point

for k = 2:npoints
    c_new = gamma_rev(:,k - 1)'; % last point
    d_new = gamma_rev(:,k - 1)'; % last point
    V_old = V_new;
    W_old = W_new;
    for i = 1:n
        c_new = c_new - W_old{i} * step_size; % update position
    end
    gamma_rev(:,k) = c_new';
    for i = 1:m
        d_new = d_new + V_old{i} * step_size * magnification; % updated position following m-dimensional field from current point
    end
    gamma_prime_rev(:,k) = d_new';
    E_new = arrayfun(@(I) (c_new - States(I,:)) / sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
    X_new = arrayfun(@(I) sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates

    W_new = cell(n,1);
    [W_new{:}] = deal(zeros(1,n));
    for i = 1:n % define ith velocity component
        for j = 1:n
            if i < j
                W_new{i} = W_new{i} + (-1) ^ (i + j - 1) * X_new{j} * E_new{i}; % velocity along ith basis vector due to jth distance
            elseif i > j
                W_new{i} = W_new{i} + (-1) ^ (i + j) * X_new{j} * E_new{i};
            end
        end
    end

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

gamma = [fliplr(gamma_rev(:,2:end)),gamma(:,2:end)]; % join trajectories
gamma_prime = [fliplr(gamma_prime_rev(:,2:end)),gamma_prime(:,2:end)]; % join trajectories
% gamma = fliplr(gamma_rev(:,2:end));
% gamma_prime = fliplr(gamma_prime_rev(:,2:end)); % join trajectories

cat(1,gamma,gamma_prime)
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
    plot3(gamma(1,:),gamma(2,:),c ^ 2 / 2 * ones(size(gamma,2)),'LineWidth',2,'Color','k')
    plot3(gamma_prime(1,:),gamma_prime(2,:),c ^ 2 / 2 * ones(size(gamma,2)),'LineWidth',2,'Color','r','LineStyle',':')
    quiver3(x{1},x{2},H,w{1},w{2},zeros(size(H)),'Color','k','LineWidth',1);
    axis([0,c,0,c,0,Inf])
    % daspect([1,1,1])
    
    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.ZLabel.String = '$$H$$';
    ax.ZLabel.Rotation = 0;
    set(get(gca,'XLabel'),'Interpreter','latex')
    set(get(gca,'YLabel'),'Interpreter','latex')
    set(get(gca,'ZLabel'),'Interpreter','latex')
    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
    hold off

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
    cb.Label.Interpreter = 'latex';
    cb.Label.Rotation = 0;
    set(cb,'FontSize',16)
    hold on
    plot3(gamma(1,:),gamma(2,:),c ^ 2 / 2 * ones(size(gamma,2)),'LineWidth',2,'Color','k')
    plot3(gamma_prime(1,:),gamma_prime(2,:),c ^ 2 / 2 * ones(size(gamma,2)),'LineWidth',2,'Color','r','LineStyle',':')
    scatter3(States(:,1),States(:,2),zeros(size(States,1),1),200,...
        'filled','r','o','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
    quiver3(x{1},x{2},H,w{1},w{2},zeros(size(H)),'Color','k','LineWidth',1);
    axis([0,c,0,c,0,Inf])
    % daspect([1,1,1])
    
    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.YLabel.Rotation = 0;
    ax.ZLabel.String = '$$H$$';
    ax.ZLabel.Rotation = 0;
    set(get(gca,'XLabel'),'Interpreter','latex')
    set(get(gca,'YLabel'),'Interpreter','latex')
    set(get(gca,'ZLabel'),'Interpreter','latex')
    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
    set(ax,'defaulttextinterpreter','latex')
    hold off

    if save
        exportgraphics(gcf,sprintf('%s/Distorted trajectory n = %d, start = (%0.6f,%0.6f), state = (%0.6f,%0.6f), magnification = %d.png',pwd,n,c0(1),c0(2),States(m,1),States(m,2),magnification),'Resolution',300)
    end
end

