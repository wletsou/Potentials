% Plots ellipses for a range of c values in a fixed grid

% PotentialField5(2,'steps',10,'npoints',500,'step_size',5e-3,'save',true)

function PotentialField5(varargin)

p = inputParser;
p.CaseSensitive = 1;

addRequired(p,'n',@(X) X >= 2)
addParameter(p,'c_max',sqrt(2),@(X) X > 0); % maximum distance
addParameter(p,'steps',10, @(X) X >= 1); % number of divisions of c
addParameter(p,'extraStates',[]); % matrix of extra reference states
addParameter(p,'startPoint',[]); % number of divisions of c
addParameter(p,'npoints',1000); % number of integration points
addParameter(p,'step_size',1e-3); % integration step size
addParameter(p,'tol',1e-3); % ellipse solution tolerance (absolute)
addParameter(p,'save',false); % whether to save image or not

parse(p,varargin{:})

n = p.Results.n; % dimension of space, initial number of reference states
c_max = p.Results.c_max; % maximum distance between reference states
steps = p.Results.steps; % number of divisions of c
npoints = p.Results.npoints; % number of integration points
step_size = p.Results.step_size; % number of integration points
tol = p.Results.tol; % ellipse solution tolerance (absolute)
save = p.Results.save;

dc = c_max / steps; 

States = eye(n) * c_max / sqrt(2); % location of fixed reference states
if ~isempty(p.Results.extraStates) && size(p.Results.extraStates,2) == n
    States = cat(1,States,p.Results.extraStates);
    m = size(States,1); % new number of reference states
else
    m = n;
end

char = strlength(num2str(m)); % number of characters in the number m of states
str = '[x1'; % create m-dimensional mesh, for vector fields
for i = 2:n
    fmt = sprintf('%%s,x%%0.%sd',num2str(char));
    str = sprintf(fmt,str,i);
end
str = sprintf('%s] = ndgrid(0:dc:c_max);',str);
eval(str);

x = cell(m,1);
[x{:}] = deal(zeros(size(x1)));
for i = 1:n % transfer X's into a cell array
    eval(sprintf('x{%d} = x%d;',i,i))
end

str = '[y1'; % create n-dimensional mesh, for ellipse points
for i = 2:n
    fmt = sprintf('%%s,y%%0.%sd',num2str(char));
    str = sprintf(fmt,str,i);
end
str = sprintf('%s] = ndgrid(0:(c_max * tol * 10):c_max);',str);
eval(str);

y = cell(n,1);
[y{:}] = deal(zeros(size(y1)));
for i = 1:n % transfer X's into a cell array
    eval(sprintf('y{%d} = y%d;',i,i))
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

H = sum(cat(n + 1,X{:}) .^ 2,n + 1) / 2; %energy function

% integral curve
trajectories = cell(steps + 1,1);
ellipses = cell(steps + 1,1);
r = 1;
for c = (0.5 * c_max):(dc * 0.5):c_max
% for c = 0:dc:c_max
    a = c_max; % ellipse major axis length
    b = c ^ 2 / c_max; % ellipse minor axis length
    eccentricity = sqrt(1 - b ^ 2 / a ^ 2);
    d = c_max * eccentricity / 2; % half the interfocal distance (n = 2), or the distance of the focus from the center

    foci = d * (eye(n) * c_max / sqrt(2) - ones(n) * c_max / 2 / sqrt(n)) ...
        / (c_max * sqrt((sqrt(2 * n) - 1) ^ 2 + (n - 1)) / 2 / sqrt(n)) ...
        + ones(n) * c_max / 2 / sqrt(n); % foci a distance d from the center along unit vectors connecting vertices with the center

    Y = cell(n,1);
    [Y{:}] = deal(zeros(size(y1)));
    f = cell(n); % basis vectors for distance between vertices and center

    for i = 1:n
        for j = 1:n
            f{i,j} = y{j} - foci(i,j); % jth component of distance from ith focus 
        end
        Y{i} = sqrt(sum(cat(n + 1,f{i,:}) .^ 2,n + 1)); % distance to the ith focus
    end

    % tol = 0.001; % tolerance for solving Y1 + ... + Yn = (n - 1) * c_max

    for i = 1:n
        coords_i = y{i}; % array of all y_i coordinates
        ellipses{r} = cat(2,ellipses{r},coords_i(abs(sum(cat(n + 1,Y{:}),n + 1) - (n - 1) * c_max) < tol)); % y_i coordinates at which Y1 + ... + Yn = c_max
    end
    % ellipses{r}

    c0 = repmat(c_max / sqrt(n) / 2,[1,n]) + repmat(c / sqrt(n) / 2,[1,n]);
    fprintf('Starting point is: [')
    fprintf('%f, ',c0(1:end - 1))
    fprintf('%f]\n,',c0(end))

    E_new = arrayfun(@(I) (c0 - States(I,:)) / sqrt(sum((c0 - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
    X_new = arrayfun(@(I) sqrt(sum((c0 - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
    
    % forward trajectory
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
    
    for k = 2:npoints
        c_new = gamma(:,k - 1)'; % last point
        V_old = V_new;
        for i = 1:m
            c_new = c_new + V_old{i} * step_size; % update position
        end
        gamma(:,k) = c_new';
        E_new = arrayfun(@(I) (c_new - States(I,:)) / sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
        X_new = arrayfun(@(I) sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
    
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
    
    for k = 2:npoints
        c_new = gamma_rev(:,k - 1)'; % last point
        V_old = V_new;
        for i = 1:m
            c_new = c_new - V_old{i} * step_size; % update position
        end
        gamma_rev(:,k) = c_new';
        E_new = arrayfun(@(I) (c_new - States(I,:)) / sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
        X_new = arrayfun(@(I) sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
    
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
    
    
    trajectories{r} = [fliplr(gamma_rev(:,2:end)),gamma]; % join left and right trajectories
    gamma = trajectories{r};
    distance = sqrt(sum((trajectories{r} - repmat(foci(:,1),[1,size(trajectories{r},2)])) .^ 2,1)) +...
        sqrt(sum((trajectories{r} - repmat(foci(:,2),[1,size(trajectories{r},2)])) .^ 2,1)); % check distance from foci (n = 2)
    % Area = cumsum(sqrt(sum((trajectories{r} - repmat(States(:,1),[1,size(trajectories{r},2)])) .^ 2,1)) + ...
    %     sqrt(sum((trajectories{r} - repmat(States(:,2),[1,size(trajectories{r},2)])) .^ 2,1))) * step_size; % area integral?
    % gamma(:,abs(gamma(1,:) - 0) < 0.002)
    distance(abs(gamma(1,:) - 0) < 0.002)

    pi * c ^ 2 / 2;
    r = r + 1;
end

if n == 2
    if ishandle(1)
        set(0,'CurrentFigure',1)
        cla
    else
        figure(1)
    end
    
    
    hold on
    for i = 1:(steps + 1)
        gamma = trajectories{i};
        ellipse = ellipses{i};
        ellipse = ellipse(ellipse(:,1) + ellipse(:,2) >= c_max / sqrt(n),:); % only get values above diagonal
        c = c_max * 0.5 + 0.5 * dc * (i - 1);
        % c = dc * (i - 1);
        plot(ellipse(:,1),ellipse(:,2),...
            'Marker','o','MarkerFaceColor',[1,0,0] * (i) / (steps + 1),'LineStyle','none','Color',[1,0,0] * (i) / (steps + 1),...
            'DisplayName',sprintf("c = %0.2f",c))
        hleg = legend('show');
        hleg.Location = 'southwest';
        plot(gamma(1,:),gamma(2,:),...
            'LineWidth',2,'LineStyle','-','Color',[1,0,0] * (i) / (steps + 1))
        hleg = legend('show');
        hleg.String(end) = [];
        hleg.Location = 'southwest';
        
        % quiver3(x{1},x{2},-H,v{1},v{2},zeros(size(H)),'Color','k','LineWidth',1);
        % plot(y1,y2,'Marker','.') % plot gridpoints
        
    end
    surf(x{1},x{2},H,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.3);
    hleg = legend('show');
    hleg.String(end) = [];
    hleg.Location = 'southwest';
    view(0,90)
    cb = colorbar;
    set(cb,'Limits',[0,max(H(:))])
    cb.Label.String = '$$H$$';
    cb.Label.Interpreter = 'latex';
    cb.Label.Rotation = 0;
    axis([0,c_max,0,c_max,0,Inf])
        % daspect([1,1,1])
    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.ZLabel.String = '$$H$$';
    zp = get(get(gca,'ZLabel'),'Position');
    zp(2) = 3 * zp(2);
    set(get(gca,'ZLabel'),'Position',zp)
    ax.ZLabel.Rotation = 0;
    ax.YLabel.Rotation = 0;

    clim([0,c_max])

    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
    set(get(gca,'XLabel'),'Interpreter','latex')
    set(get(gca,'YLabel'),'Interpreter','latex')
    set(get(gca,'ZLabel'),'Interpreter','latex')
    hold off

    if save
        exportgraphics(gcf,sprintf('%s/Deformed trajectory series n = %d.png',pwd,n),'Resolution',300)
    end

end

