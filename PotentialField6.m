% Plot a range of (disorted) two-dimensional trajectories and their
% undistorted counterparts

% PotentialField6(2,'steps',10,'npoints',500,'step_size',5e-3,'save',true)

function PotentialField6(varargin)

p = inputParser;
p.CaseSensitive = 1;

addRequired(p,'n',@(X) X == 2)
addParameter(p,'c_max',sqrt(2),@(X) X > 0); % maximum distance
addParameter(p,'steps',10, @(X) X >= 1); % number of divisions of c
addParameter(p,'extraStates',[]); % matrix of extra reference states
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
str = sprintf('%s] = ndgrid(0:(c_max / 20):c_max);',str);
eval(str);

x = cell(m,1);
[x{:}] = deal(zeros(size(x1)));
for i = 1:n % transfer X's into a cell array
    eval(sprintf('x{%d} = x%d;',i,i))
end

str = '[y1'; % create n-dimensional mesh, for alternative reference states points
for i = 2:n
    fmt = sprintf('%%s,y%%0.%sd',num2str(char));
    str = sprintf(fmt,str,i);
end
str = sprintf('%s] = ndgrid(0:dc:c_max);',str);
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

V = cell(m,1); % relative frame veLocities
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

v = cell(n,1); % fixed frame veLocities
[v{:}] = deal(zeros(size(x1)));
for i = 1:n
    for j = 1:m
        v{i} = v{i} + V{j} .* e{j,i}; % ith fixed frame velocity due to jth relative frame velocity
    end
end

Pairs = nchoosek(1:m,2); % Pairs of reference states
H_surf = sum(cat(n + 1,X{:}) .^ 2,n + 1) / 2; % energy surface

% innitialize integral curves and energies for each value of c
trajectories = cell(steps + 1,2); % relative to reference states, full and partial curves
energies = cell(steps + 1,2);

trajectories_alt = cell(steps + 1,1); % relative to alternative reference states
energies_alt = cell(steps + 1,1);

r = 1;
for c = 0:dc:c_max
    c_star = sqrt(c * c_max); % effective diamter

    Loci = (c_star / 2) * (eye(n) * c_max / sqrt(2) - ones(n) * c_max / n / sqrt(2)) ...
        / (c_max * sqrt((n - 1) / 2 / n)) ...
        + ones(n) * c_max / n / sqrt(2); % Loci a distance c_star / 2 from the center along unit vectors connecting vertices with the center

    Y = cell(n,1);
    [Y{:}] = deal(zeros(size(y1)));
    f = cell(n); % basis vectors for distance between vertices and center

    for i = 1:n
        for j = 1:n
            f{i,j} = y{j} - Loci(i,j); % jth component of distance from ith locus 
        end
        Y{i} = sqrt(sum(cat(n + 1,f{i,:}) .^ 2,n + 1)); % distance to the ith locus
    end

    c0 = repmat(c_max / sqrt(2) / n,[1,n]) + repmat(c / sqrt(2) / n,[1,n]);
    fprintf('Starting point is: [')
    fprintf('%f, ',c0(1:end - 1))
    fprintf('%f]\n,',c0(end))

    % forward trajectory
    E_new = arrayfun(@(I) (c0 - States(I,:)) / sqrt(sum((c0 - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
    X_new = arrayfun(@(I) sqrt(sum((c0 - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
    Phi_new = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(m,2),'UniformOutput',false);
    Phi_new = sum(cat(1,Phi_new{:}),1); % potential function
    H_new = arrayfun(@(I) X_new{I} .^ 2 / 2,1:m,'UniformOutput',false);
    H_new = sum(cat(1,H_new{:}),1); % energy function

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
    Phi = zeros(1,npoints);
    Phi(1) = Phi_new;
    H = zeros(1,npoints); % initialize vector of computed energies
    H(1) = H_new;

    for k = 2:npoints
        c_new = gamma(:,k - 1)'; % last point
        V_old = V_new;
        for i = 1:m
            c_new = c_new + V_old{i} * step_size; % update position
        end
        gamma(:,k) = c_new';
        E_new = arrayfun(@(I) (c_new - States(I,:)) / sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
        X_new = arrayfun(@(I) sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
        
        % update functions along trajectory
        Phi_new = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(m,2),'UniformOutput',false);
        Phi_new = sum(cat(1,Phi_new{:}),1); % potential function
        Phi(k) = Phi_new;
        H_new = arrayfun(@(I) X_new{I} .^ 2 / 2,1:m,'UniformOutput',false);
        H_new = sum(cat(1,H_new{:}),1); % energy function
        H(k) = H_new;

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

    % cumulative distance, trajectory points within a distance pi * c_star / 4 of center
    gamma_sub = gamma(:,[true,cumsum(sqrt(sum(diff(gamma,1,2) .^ 2,1))) <= pi * c_star / 4]);
    H_sub = H(:,[true,cumsum(sqrt(sum(diff(gamma,1,2) .^ 2,1))) <= pi * c_star / 4]);

    % reverse trajectory
    E_new = arrayfun(@(I) (c0 - States(I,:)) / sqrt(sum((c0 - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
    X_new = arrayfun(@(I) sqrt(sum((c0 - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
    Phi_new_rev = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(m,2),'UniformOutput',false);
    Phi_new_rev = sum(cat(1,Phi_new_rev{:}),1); % potential function
    H_new_rev = arrayfun(@(I) X_new{I} .^ 2 / 2,1:m,'UniformOutput',false);
    H_new_rev = sum(cat(1,H_new_rev{:}),1); % energy function

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
    Phi_rev = zeros(1,npoints);
    Phi_rev(1) = Phi_new_rev;
    H_rev = zeros(1,npoints); % initialize vector of computed energies
    H_rev(1) = H_new_rev;

    for k = 2:npoints
        c_new = gamma_rev(:,k - 1)'; % last point
        V_old = V_new;
        for i = 1:m
            c_new = c_new - V_old{i} * step_size; % update position
        end
        gamma_rev(:,k) = c_new';

        % update position basis vectors
        E_new = arrayfun(@(I) (c_new - States(I,:)) / sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
        X_new = arrayfun(@(I) sqrt(sum((c_new - States(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
        
        % update functions along trajectory
        Phi_new_rev = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(m,2),'UniformOutput',false);
        Phi_new_rev = sum(cat(1,Phi_new_rev{:}),1); % potential function
        Phi_rev(k) = Phi_new_rev;
        H_new_rev = arrayfun(@(I) X_new{I} .^ 2 / 2,1:m,'UniformOutput',false);
        H_new_rev = sum(cat(1,H_new_rev{:}),1); % energy function
        H_rev(k) = H_new_rev;
        
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

    % cumulative distance, trajectory points within a distance pi * c_star / 4 of center
    gamma_rev_sub = gamma_rev(:,[true,cumsum(sqrt(sum(diff(gamma_rev,1,2) .^ 2,1))) <= pi * c_star / 4]);
    H_rev_sub = H_rev(:,[true,cumsum(sqrt(sum(diff(gamma_rev,1,2) .^ 2,1))) <= pi * c_star / 4]); 

    trajectories{r,1} = [fliplr(gamma_rev(:,2:end)),gamma]; % join left and right trajectories
    trajectories{r,2} = [fliplr(gamma_rev_sub(:,2:end)),gamma_sub]; % join left and right trajectories of partial curves
    energies{r,1} = [fliplr(H_rev(:,2:end)),H]; % join left and right trajectories
    energies{r,2} = [fliplr(H_rev_sub(:,2:end)),H_sub]; % join left and right trajectories of partial curves

    % alternative trajectories with separation distance c_star
    c0 = repmat(c_max / sqrt(2) / n,[1,n]) + repmat(c_star / sqrt(2) / n,[1,n]);
    fprintf('Alternative starting point is: [')
    fprintf('%f, ',c0(1:end - 1))
    fprintf('%f]\n,',c0(end))

    % forward trajectory
    E_new = arrayfun(@(I) (c0 - Loci(I,:)) / sqrt(sum((c0 - Loci(I,:)) .^ 2)),1:n,'UniformOutput',false); % unit vectors toward alternative reference states
    X_new = arrayfun(@(I) sqrt(sum((c0 - Loci(I,:)) .^ 2)),1:n,'UniformOutput',false); % relative frame coordinates
    Phi_new = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(n,2),'UniformOutput',false);
    Phi_new = sum(cat(1,Phi_new{:}),1); % potential function
    H_new = arrayfun(@(I) X_new{I} .^ 2 / 2,1:n,'UniformOutput',false);
    H_new = sum(cat(1,H_new{:}),1); % energy function

    V_new = cell(n,1);
    [V_new{:}] = deal(zeros(1,n));
    for i = 1:n % define ith (alternative) velocity component
        for j = 1:n
            if i < j
                V_new{i} = V_new{i} + (-1) ^ (i + j - 1) * X_new{j} * E_new{i}; % velocity along ith basis vector due to jth distance
            elseif i > j
                V_new{i} = V_new{i} + (-1) ^ (i + j) * X_new{j} * E_new{i};
            end
        end
    end
    
    gamma = zeros(n,npoints); % initialize matrix of sample points
    gamma(:,1) = c0'; % first point
    Phi = zeros(1,npoints);
    Phi(1) = Phi_new;
    H = zeros(1,npoints); % initialize vector of computed energies
    H(1) = H_new;

    for k = 2:npoints
        c_new = gamma(:,k - 1)'; % last point
        V_old = V_new;
        for i = 1:n
            c_new = c_new + V_old{i} * step_size; % update position
        end
        gamma(:,k) = c_new';
        E_new = arrayfun(@(I) (c_new - Loci(I,:)) / sqrt(sum((c_new - Loci(I,:)) .^ 2)),1:n,'UniformOutput',false); % unit vectors toward alternative reference states
        X_new = arrayfun(@(I) sqrt(sum((c_new - Loci(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
        
        % update functions along trajectory
        Phi_new = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(n,2),'UniformOutput',false);
        Phi_new = sum(cat(1,Phi_new{:}),1); % potential function
        Phi(k) = Phi_new;
        H_new = arrayfun(@(I) X_new{I} .^ 2 / 2,1:n,'UniformOutput',false);
        H_new = sum(cat(1,H_new{:}),1); % energy function
        H(k) = H_new;

        V_new = cell(n,1);
        [V_new{:}] = deal(zeros(1,n));
        for i = 1:n % define ith velocity component
            for j = 1:n
                if i < j
                    V_new{i} = V_new{i} + (-1) ^ (i + j - 1) * X_new{j} * E_new{i}; % velocity along ith basis vector due to jth distance
                elseif i > j
                    V_new{i} = V_new{i} + (-1) ^ (i + j) * X_new{j} * E_new{i};
                end
            end
        end
    end
    
    % reverse trajectory
    E_new = arrayfun(@(I) (c0 - Loci(I,:)) / sqrt(sum((c0 - Loci(I,:)) .^ 2)),1:n,'UniformOutput',false); % unit vectors toward alternative reference states
    X_new = arrayfun(@(I) sqrt(sum((c0 - Loci(I,:)) .^ 2)),1:n,'UniformOutput',false); % relative frame coordinates
    Phi_new_rev = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(n,2),'UniformOutput',false);
    Phi_new_rev = sum(cat(1,Phi_new_rev{:}),1); % potential function
    H_new_rev = arrayfun(@(I) X_new{I} .^ 2 / 2,1:m,'UniformOutput',false);
    H_new_rev = sum(cat(1,H_new_rev{:}),1); % energy function

    V_new = cell(n,1);
    [V_new{:}] = deal(zeros(1,n));
    for i = 1:n % define ith velocity component
        for j = 1:n
            if i < j
                V_new{i} = V_new{i} + (-1) ^ (i + j - 1) * X_new{j} * E_new{i}; % velocity along ith basis vector due to jth distance
            elseif i > j
                V_new{i} = V_new{i} + (-1) ^ (i + j) * X_new{j} * E_new{i};
            end
        end
    end
    
    gamma_rev = zeros(n,npoints); % initialize matrix of sample points
    gamma_rev(:,1) = c0'; % first point
    Phi_rev = zeros(1,npoints);
    Phi_rev(1) = Phi_new_rev;
    H_rev = zeros(1,npoints); % initialize vector of computed energies
    H_rev(1) = H_new_rev;

    for k = 2:npoints
        c_new = gamma_rev(:,k - 1)'; % last point
        V_old = V_new;
        for i = 1:n
            c_new = c_new - V_old{i} * step_size; % update position
        end
        gamma_rev(:,k) = c_new';

        % update position basis vectors
        E_new = arrayfun(@(I) (c_new - Loci(I,:)) / sqrt(sum((c_new - Loci(I,:)) .^ 2)),1:n,'UniformOutput',false); % unit vectors toward alternative reference states
        X_new = arrayfun(@(I) sqrt(sum((c_new - Loci(I,:)) .^ 2)),1:n,'UniformOutput',false); % relative frame coordinates
        
        % update functions along trajectory
        Phi_new_rev = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(n,2),'UniformOutput',false);
        Phi_new_rev = sum(cat(1,Phi_new_rev{:}),1); % potential function
        Phi_rev(k) = Phi_new_rev;
        H_new_rev = arrayfun(@(I) X_new{I} .^ 2 / 2,1:n,'UniformOutput',false);
        H_new_rev = sum(cat(1,H_new_rev{:}),1); % energy function
        H_rev(k) = H_new_rev;
        
        V_new = cell(n,1);
        [V_new{:}] = deal(zeros(1,n));
        for i = 1:n % define ith velocity component
            for j = 1:n
                if i < j
                    V_new{i} = V_new{i} + (-1) ^ (i + j - 1) * X_new{j} * E_new{i}; % velocity along ith basis vector due to jth distance
                elseif i > j
                    V_new{i} = V_new{i} + (-1) ^ (i + j) * X_new{j} * E_new{i};
                end
            end
        end
    end
    trajectories_alt{r} = [fliplr(gamma_rev(:,2:end)),gamma]; % join left and right trajectories
    energies_alt{r} = [fliplr(H_rev(:,2:end)),H]; % join left and right trajectories

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
        gamma = trajectories{i,1}; % full trajectory
        gamma_sub = trajectories{i,2}; % partial trajectory
        gamma_alt = trajectories_alt{i}; % corresponding alternative trajectory
        h = energies{i,1}; % energy computed along full trajectory
        h_sub = energies{i,2}; % energy computed along full trajectory
        h_alt = energies_alt{i}; % energy computed along corresponding alternative trajectory

        c = dc * (i - 1);
        c_star = sqrt(c_max * c);

        plot3(gamma(1,:),gamma(2,:),h,...
            'LineWidth',2,'LineStyle','-','Color',[1,0,0] * (i) / (steps + 1),...
            'DisplayName',sprintf("c = %0.2f, c* = %0.2f",c,c_star))
        hleg = legend('show');
        hleg.Location = 'southwest';
        
        plot3(gamma_sub(1,:),gamma_sub(2,:),h_sub,...
            'LineWidth',4,'LineStyle','-','Color',[1,0,0] * (i) / (steps + 1))
        hleg = legend('show');
        hleg.String(end) = [];

        plot3(gamma_alt(1,:),gamma_alt(2,:),h_alt + (c_max ^ 2 - c_star ^ 2) / 4,...
            'LineWidth',2,'LineStyle','--','Color',[1,0,0] * (i) / (steps + 1))
        hleg = legend('show');
        hleg.String(end) = [];        
    end
    % uncomment to plot vector field relative to fixed reference states
    % quiver3(x{1},x{2},H_surf,v{1},v{2},zeros(size(H)),'Color','k','LineWidth',1);
    % hleg = legend('show');
    % hleg.String(end) = [];

    surf(x{1},x{2},H_surf,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.3); % plot energy surface
    hleg = legend('show');
    hleg.String(end) = [];
    hleg.Location = 'southeast';

    grid on

    view(20,30)

    % cb = colorbar;
    % set(cb,'Limits',[0,max(H_surf(:))])
    % cb.Label.String = '$$H$$';
    % cb.Label.Interpreter = 'latex';
    % cb.Label.Rotation = 0;
    % set(cb,'FontSize',16)
    clim([c_max ^ 2 / 4,c_max ^ 2 / 2])

    axis([0,c_max,0,c_max,0.95 * c_max ^ 2 / 4,1.05 * c_max ^ 2 / 2])
    % daspect([1,1,1])
    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.ZLabel.String = '$$H$$';
    zp = get(get(gca,'ZLabel'),'Position');
    zp(2) = 5 * zp(2);
    set(get(gca,'ZLabel'),'Position',zp)
    ax.ZLabel.Rotation = 0;
    % ax.YLabel.Rotation = 0;

    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
    set(get(gca,'XLabel'),'Interpreter','latex')
    set(get(gca,'YLabel'),'Interpreter','latex')
    set(get(gca,'ZLabel'),'Interpreter','latex')
    hold off

    if save
        exportgraphics(gcf,sprintf('%s/Deformed trajectory series n = %d.png',pwd,n),'Resolution',300)
    end


    if ishandle(2)
        set(0,'CurrentFigure',2)
        cla
    else
        figure(2)
    end
    
    hold on
    for i = 1:(steps + 1)
        gamma = trajectories{i,1}; % full trajectory
        gamma_sub = trajectories{i,2}; % partial trajectory
        gamma_alt = trajectories_alt{i}; % corresponding alternative trajectory
        h = energies{i,1}; % energy computed along full trajectory
        h_sub = energies{i,2}; % energy computed along full trajectory
        h_alt = energies_alt{i}; % energy computed along corresponding alternative trajectory

        c = dc * (i - 1);
        c_star = sqrt(c_max * c);

        plot3(gamma(1,:),gamma(2,:),h,...
            'LineWidth',2,'LineStyle','-','Color',[1,0,0] * (i) / (steps + 1),...
            'DisplayName',sprintf("c = %0.2f, c* = %0.2f",c,c_star))
        hleg = legend('show');
        hleg.Location = 'southwest';

        plot3(gamma_sub(1,:),gamma_sub(2,:),h_sub,...
            'LineWidth',4,'LineStyle','-','Color',[1,0,0] * (i) / (steps + 1))
        hleg = legend('show');
        hleg.String(end) = [];
        
        plot3(gamma_alt(1,:),gamma_alt(2,:),h_alt + (c_max ^ 2 - c_star ^ 2) / 4,...
            'LineWidth',2,'LineStyle','--','Color',[1,0,0] * (i) / (steps + 1))
        hleg = legend('show');
        hleg.String(end) = [];
    end
    % uncomment to plot vector field relative to fixed reference states
    % quiver3(x{1},x{2},H_surf,v{1},v{2},zeros(size(H)),'Color','k','LineWidth',1);
    % hleg = legend('show');
    % hleg.String(end) = [];

    surf(x{1},x{2},H_surf,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.3); % plot energy surface
    hleg = legend('show');
    hleg.String(end) = [];
    hleg.Location = 'southwest';

    grid on

    view(0,90)

    cb = colorbar;
    set(cb,'Limits',[0,max(H_surf(:))])
    cb.Label.String = '$$H$$';
    cb.Label.Interpreter = 'latex';
    cb.Label.Rotation = 0;
    set(cb,'FontSize',16)
    clim([c_max ^ 2 / 4,c_max ^ 2 / 2])

    axis([0,c_max,0,c_max,0,Inf])
    % daspect([1,1,1])
    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.ZLabel.String = '$$H$$';
    zp = get(get(gca,'ZLabel'),'Position');
    zp(2) = 5 * zp(2);
    set(get(gca,'ZLabel'),'Position',zp)
    ax.ZLabel.Rotation = 0;
    ax.YLabel.Rotation = 0;

    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
    set(get(gca,'XLabel'),'Interpreter','latex')
    set(get(gca,'YLabel'),'Interpreter','latex')
    set(get(gca,'ZLabel'),'Interpreter','latex')
    hold off

    if save
        exportgraphics(gcf,sprintf('%s/Deformed trajectory series n = %d, top-down view.png',pwd,n),'Resolution',300)
    end
end

