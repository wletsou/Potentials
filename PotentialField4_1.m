% Plots trajectory (no vector field) for a range of c values in a fixed
% grid, shifting origin

% PotentialField4_1(2,'steps',10,'npoints',200,'step_size',5e-3,'save',true)

function PotentialField4_1(varargin)

p = inputParser;
p.CaseSensitive = 1;

addRequired(p,'n',@(X) X >= 2)
addParameter(p,'c_max',sqrt(2),@(X) X > 0); % maximum distance
addParameter(p,'steps',10, @(X) X >= 1); % number of divisions of c
addParameter(p,'extraStates',[]); % matrix of extra reference states
addParameter(p,'startPoint',[]); % number of divisions of c
addParameter(p,'npoints',1000); % number of integration points
addParameter(p,'step_size',1e-3); % integration step size
addParameter(p,'save',false); % whether to save image or not

parse(p,varargin{:})

n = p.Results.n; % dimension of space, initial number of reference states
c_max = p.Results.c_max; % maximum distance between reference states
steps = p.Results.steps; % number of divisions of c
npoints = p.Results.npoints; % number of integration points
step_size = p.Results.step_size; % number of integration points
save = p.Results.save;

dc = c_max / steps; 

trajectories = cell(steps + 0,1);
energies = cell(steps + 0,1);
r = 1;
for c = 0:dc:c_max

    States = (eye(n) * (c_max + c) + (ones(n) - eye(n)) * (c_max - c)) / sqrt(2) / 2 - ones(n) * c / 2 / sqrt(2); % shift endpoints so that center of curve is always in the same place
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
    
    H = sum(cat(n + 1,X{:}) .^ 2,n + 1) / 2; %energy function

    % integral curve

    c0 = repmat(c_max / sqrt(n) / 2,[1,n]); % always start midway between the reference states
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
    energy = zeros(1,npoints);
    energy(1) = sum(cat(n + 1,X_new{:}) .^ 2,n + 1) / 2;
    
    for k = 2:npoints
        c_new = gamma(:,k - 1)'; % last point
        V_old = V_new;
        for i = 1:m
            c_new = c_new + V_old{i} * step_size; % update position
        end
        gamma(:,k) = c_new';
        E_new = arrayfun(@(I) (c_new - States(I,:)) / sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
        X_new = arrayfun(@(I) sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
        
        energy(k) = sum(cat(n + 1,X_new{:}) .^ 2,n + 1) / 2;

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
    energy_rev = zeros(1,npoints);
    energy_rev(1) = sum(cat(n + 1,X_new{:}) .^ 2,n + 1) / 2;

    
    for k = 2:npoints
        c_new = gamma_rev(:,k - 1)'; % last point
        V_old = V_new;
        for i = 1:m
            c_new = c_new - V_old{i} * step_size; % update position
        end
        gamma_rev(:,k) = c_new';
        E_new = arrayfun(@(I) (c_new - States(I,:)) / sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
        X_new = arrayfun(@(I) sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
    
        energy_rev(k) = sum(cat(n + 1,X_new{:}) .^ 2,n + 1) / 2;

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
    energies{r} = [fliplr(energy_rev(2:end)),energy];
    r = r + 1;
end

if n == 2
    if ishandle(1)
        set(0,'CurrentFigure',1)
        cla
    else
        figure(1)
    end
    
    combined_trajectories = cat(1,trajectories{:});
    combined_energies = cat(1,energies{:});

    xdat = combined_trajectories(1:2:end,:);
    ydat = combined_trajectories(2:2:end,:);
    surf(xdat(:,1:10:end),ydat(:,1:10:end),-combined_energies(1:end,1:10:end),'EdgeColor','none','FaceColor','interp','FaceAlpha',0.3);
    view(30,30)
    hold on
    for i = 1:(steps+1)
        
        gamma = trajectories{i};
        energy = energies{i};
        c = dc * (i - 1);
        States = (eye(n) * (c_max + c) + (ones(n) - eye(n)) * (c_max - c)) / sqrt(2) / 2 - ones(n) * c / 2 / sqrt(2);
        States = cat(1,States,ones(n) * c_max / 2 / sqrt(2));
        plot3(gamma(1,:),gamma(2,:),-energy,...
            'LineWidth',2,'LineStyle','-','Color',[1,0,1] * (steps + 1 - i + 1) / (steps + 1))
        plot3(gamma(1,:),gamma(2,:),-energy,...
            'LineWidth',2,'LineStyle','-','Color',[1,0,1] * (steps + 1 - i + 1) / (steps + 1))
        % quiver3(x{1},x{2},-H,v{1},v{2},zeros(size(H)),'Color','k','LineWidth',1);
        scatter3(States(1,1),States(1,2),-c ^ 2 / 2 * ones(size(States,1),1),200,...
            'filled','b','o','MarkerEdgeColor','none','MarkerFaceAlpha',(i) / (steps + 1))
        scatter3(States(2,1),States(2,2),-c ^ 2 / 2 * ones(size(States,1),1),200,...
            'filled','r','o','MarkerEdgeColor','none','MarkerFaceAlpha',(i) / (steps + 1))
        scatter3(States(3,1),States(3,2),-c ^ 2 / 2 * ones(size(States,1),1),200,...
            'filled','m','o','MarkerEdgeColor','none','MarkerFaceAlpha',(i) / (steps + 1))
        
    end
    axis([0,c,0,c,-c_max ^ 2 / 2 - 0.05,0 + 0.05] - [1,1,1,1,0,0] * c / 2 / sqrt(2))
        % daspect([1,1,1])
    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.ZLabel.String = '$$-H$$';
    zp = get(get(gca,'ZLabel'),'Position');
    zp(2) = 1.4 * zp(2);
    set(get(gca,'ZLabel'),'Position',zp)
    ax.ZLabel.Rotation = 0;

    clim([-c_max ^ 2 / 2,0])

    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
    set(get(gca,'XLabel'),'Interpreter','latex')
    set(get(gca,'YLabel'),'Interpreter','latex')
    set(get(gca,'ZLabel'),'Interpreter','latex')
    hold off

    if save
        exportgraphics(gcf,sprintf('%s/Trajectory series (shifted) n = %d.png',pwd,n),'Resolution',300)
    end

end

