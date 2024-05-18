function PotentialField11(varargin)
% Plot 2d or 3d surface with apex of sphere reduced by a factor c / c_max
% and show a vector field through apex, with the apex replacing the nth
% reference state

% PotentialField11(3,sqrt(2),'steps',100,'npoints',3000,'step_size',1e-3,'save',true)
% PotentialField11(3,1,'steps',100,'npoints',3000,'step_size',1e-3,'save',true)
% PotentialField11(3,0.5,'steps',100,'npoints',3000,'step_size',1e-3,'save',true)
% PotentialField11(3,0.2,'steps',100,'npoints',3000,'step_size',1e-3,'save',true)

rng(20200116,"twister")

p = inputParser;
p.CaseSensitive = 1;

addRequired(p,'n',@(X) X >= 2) % number of states
addRequired(p,'c',@(X) X > 0); % (scaled) distance from center of simplex containing the n states
addParameter(p,'c_max',sqrt(2),@(X) X > 0); % maximum separation distance distance
addParameter(p,'steps',10, @(X) X >= 1); % number of divisions of c
addParameter(p,'extraStates',[]); % matrix of extra reference states
addParameter(p,'npoints',1000); % number of integration points
addParameter(p,'step_size',1e-3); % integration step size
addParameter(p,'save',false); % whether to save image or not

parse(p,varargin{:})

n = p.Results.n; % dimension of space, initial number of reference states
c = p.Results.c; % scaled) distance from center of simplex containing the n states
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
str = sprintf('%s] = ndgrid(-c_max:(2 * dc):c_max);',str);
eval(str);

x = cell(m,1);
[x{:}] = deal(zeros(size(x1)));
for i = 1:n % transfer X's into a cell array
    eval(sprintf('x{%d} = x%d;',i,i))
end

% define rotation matrices for converting rectangular coordinates to
% parallel/perpedicular components
if n == 2
    R = [cos(-pi / 4) -sin(-pi / 4); sin(-pi / 4) cos(-pi / 4)]; % 2d rotation matrix
elseif n == 3
    %3d rotation matrix: rotate about x, then z to position z at (1,1,1)
    R = [cos(-pi / 4) -sin(-pi / 4) 0; sin(-pi / 4) cos(-pi / 4) 0; 0 0 1] ...
        * [1 0 0; 0 cos(-acos(1/sqrt(3))) -sin(-acos(1/sqrt(3))); 0 sin(-acos(1/sqrt(3))) cos(-acos(1/sqrt(3)))];
end
RT = R'; % metric tensor

y = cell(m,1); % transformed coordinates
[y{:}] = deal(zeros(size(x1)));

for i = 1:size(R,1)
    for j = 1:size(R,2)
        y{i} = y{i} + RT(i,j) * x{j};
    end
end

% stretch perpendicular component
y{n} = (y{n} - c_max / sqrt(2 * n)) * c_max / c + c_max / sqrt(2 * n);

% back-transform
z = cell(m,1); % un-transformed coordinates
[z{:}] = deal(zeros(size(x1)));
for i = 1:size(R,2)
    for j = 1:size(R,1)
        z{i} = z{i} + y{j} * RT(j,i);
    end
end

% compute (X ^ 2 + Y ^ 2 + ...) / 2 for each point in the un-transformed space
H = zeros(size(x1));
for i = 1:m
    for j = 1:n
        H = H + (z{j} - States(i,j)) .^ 2 / 2;
    end
end

if n == 2

    if ishandle(1)
        set(0,'CurrentFigure',1)
        cla
    else
        figure(1)
    end

    tol = 0.05;
    p1 = x{1}(abs(H - (n - 1) * c_max ^ 2 / 2) < tol); % find x, y coordinates corresponding to set value of H = c_max ^ 2 / 2
    p2 = x{2}(abs(H - (n - 1) * c_max ^ 2 / 2) < tol);
    scatter(p1,p2,'filled','b','o','MarkerEdgeColor','none','MarkerFaceAlpha',1)
    % hleg = legend('show');
    % hleg.String(end) = [];
    % hleg.Location = 'southeast';
    axis([0,1.2 * c_max / sqrt(2),0,1.2 * c_max / sqrt(2)])
    axis([-1.2 * c_max / sqrt(2),1.2 * c_max / sqrt(2),-1.2 * c_max / sqrt(2),1.2 * c_max / sqrt(2)])

elseif n == 3

    if ishandle(1)
        set(0,'CurrentFigure',1)
        cla
    else
        figure(1)
    end
  
    tol = 0.05;
    p1 = x1(abs(H - (n - 1) * c_max ^ 2 / 2) < tol); % find x, y, z coordinates corresponding to set value of H = c_max ^ 2 / 2
    p2 = x2(abs(H - (n - 1) * c_max ^ 2 / 2) < tol);
    p3 = x3(abs(H - (n - 1) * c_max ^ 2 / 2) < tol);
    scatter3(p1,p2,p3,'filled','b','o','MarkerEdgeColor','none','MarkerFaceAlpha',0.05)
    % hleg = legend('show');
    % hleg.String(end) = [];
    % hleg.Location = 'southeast';
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

v = cell(n,1); % fixed frame veLocities
[v{:}] = deal(zeros(size(x1)));
for i = 1:n
    for j = 1:m
        v{i} = v{i} + V{j} .* e{j,i}; % ith fixed frame velocity due to jth relative frame velocity
    end
end

Pairs = nchoosek(1:m,2); % Pairs of reference states

States_old = States;
States(n,:) = ones(1,n) * (c_max / sqrt(2 * n) + sqrt( (n - 1) / (2 * n) ) * c) / sqrt(n); % replace nth state with apex
c0 = States(n,:) + (rand(1,n) - 0.5) * 0.01; % shift slightly to avoid singularity

fprintf('Apex is: [')
fprintf('%f, ',States(n,1:end - 1))
fprintf('%f]\n',States(n,end))

fprintf('Starting point is: [')
fprintf('%f, ',c0(1:end - 1))
fprintf('%f]\n',c0(end))

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

    if ~any(cellfun(@(I) abs(I - 0) < 0.01,X_new(1:n-1)) == 1) % only update if not at one of the endpoints
        V_old = V_new;
        for i = 1:m
            c_new = c_new + V_old{i} * step_size; % update position
        end
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
    if ~any(cellfun(@(I) abs(I - 0) < 0.01,X_new(1:n-1)) == 1) % only update if not at one of the endpoints
        V_old = V_new;
        for i = 1:m
            c_new = c_new - V_old{i} * step_size; % update position
        end
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

trajectories{1} = [fliplr(gamma_rev(:,2:end)),gamma]; % join left and right trajectories
energies{1} = [fliplr(H_rev(:,2:end)),H]; % join left and right trajectories

% trajectory through original three states
c0 = States_old(n,:) + (rand(1,n) - 0.5) * 0.01; % shift slightly to avoid singularity

fprintf('\nStarting point is: [')
fprintf('%f, ',c0(1:end - 1))
fprintf('%f]\n',c0(end))

% forward trajectory
E_new = arrayfun(@(I) (c0 - States_old(I,:)) / sqrt(sum((c0 - States_old(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
X_new = arrayfun(@(I) sqrt(sum((c0 - States_old(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
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
    E_new = arrayfun(@(I) (c_new - States_old(I,:)) / sqrt(sum((c_new - States_old(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
    X_new = arrayfun(@(I) sqrt(sum((c_new - States_old(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
   
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

% reverse trajectory
E_new = arrayfun(@(I) (c0 - States_old(I,:)) / sqrt(sum((c0 - States_old(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
X_new = arrayfun(@(I) sqrt(sum((c0 - States_old(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
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
    E_new = arrayfun(@(I) (c_new - States_old(I,:)) / sqrt(sum((c_new - States_old(I,:)) .^ 2)),1:m,'UniformOutput',false); % unit vectors toward reference states
    X_new = arrayfun(@(I) sqrt(sum((c_new - States_old(I,:)) .^ 2)),1:m,'UniformOutput',false); % relative frame coordinates
    
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

trajectories{2} = [fliplr(gamma_rev(:,2:end)),gamma]; % join left and right trajectories
energies{2} = [fliplr(H_rev(:,2:end)),H]; % join left and right trajectories

if n == 3
    
    hold on
    gamma = trajectories{1}; % trajectory through new state
    plot3(gamma(1,:),gamma(2,:),gamma(3,:),...
        'LineWidth',2,'LineStyle','-','Color',[1,0,0],...
        'DisplayName',sprintf("c = %0.2f",c))
    %hleg = legend('show');

    gamma = trajectories{2}; % trajectory through original states
    plot3(gamma(1,:),gamma(2,:),gamma(3,:),...
        'LineWidth',2,'LineStyle','-','Color',[0,0,1],...
        'DisplayName',sprintf("c = %0.2f",c))
    %hleg = legend('show');

    
    % scatter3(States_old(:,1),States_old(:,2),States_old(:,3),200,...
    % 'filled','o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
    
    States_old = [States_old; States_old(1,:)]; % make a closed loop
    plot3(States_old(:,1),States_old(:,2),States_old(:,3),'-','LineWidth',2,'Color','k')
    %hleg = legend('show');
    scatter3(States(3,1),States(3,2),States(3,3),200,...
    'filled','o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor','none','MarkerFaceAlpha',0.5) % new state 

    %hleg = legend('show');
    %hleg.String(end) = []; 

    view(135,45)
    view(110,25)
    view(65,15)
    axis square
    axis([0,1.2 * c_max / sqrt(2),0,1.2 * c_max / sqrt(2),0,1.2 * c_max / sqrt(2)])
    axis([-0.6 * c_max / sqrt(2),1.2 * c_max / sqrt(2),...
        -0.6 * c_max / sqrt(2),1.2 * c_max / sqrt(2),...
        -0.6 * c_max / sqrt(2),1.2 * c_max / sqrt(2)])

    title(sprintf('c^* = %0.3f\t(c = %0.3f, c_{max} = %0.3f)',(c * c_max ^ (n - 1)) ^ (1 / n),c,c_max),'FontSize',18)
    
    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.ZLabel.String = '$$z$$';
    % zp = get(get(gca,'ZLabel'),'Position');
    % zp(2) = 5 * zp(2);
    % set(get(gca,'ZLabel'),'Position',zp)
    ax.ZLabel.Rotation = 0;
    % ax.YLabel.Rotation = 0;

    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
    set(get(gca,'XLabel'),'Interpreter','latex')
    set(get(gca,'YLabel'),'Interpreter','latex')
    set(get(gca,'ZLabel'),'Interpreter','latex')

    if save
        exportgraphics(gcf,sprintf('%s/Distorted surface n = %d c= %0.3f.png',pwd,n,c),'Resolution',300)
    end
end
