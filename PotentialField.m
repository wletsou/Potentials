% Computes n-dimensional vector fields for potential system with n
% reference states, starting in the relative frame

% PotentialField(2,'steps',20,'startPoint',[1.000860e+00, 7.523451e-04])
% PotentialField(3,'steps',10,'startPoint',[0.7930,0.2232,0.9758],'npoints',5000,'step_size',1e-3,'save',true)

function PotentialField(varargin)

p = inputParser;
p.CaseSensitive = 1;

addRequired(p,'n',@(X) X >= 2)
addParameter(p,'c',sqrt(2),@(X) X > 0); % distance between each pair of states
addParameter(p,'steps',10, @(X) X >= 1); % number of divisions of c
addParameter(p,'startPoint',[]); % number of divisions of c
addParameter(p,'npoints',1000); % number of integration points
addParameter(p,'step_size',1e-3); % integration step size
addParameter(p,'save',false); % whether to save image or not

parse(p,varargin{:})

n = p.Results.n; % number of reference states
c = p.Results.c; % distance between reference states
steps = p.Results.steps; % number of divisions of c
npoints = p.Results.npoints; % number of integration points
step_size = p.Results.step_size; % number of integration points
save = p.Results.save; % whether to save image or not

dc = c / steps; 

States = eye(n) * c / sqrt(2); % coordinates of reference states

char = strlength(num2str(n)); % number of characters in the number n of states
str = '[X1'; % create n-dimensional mesh
for i = 2:n
    fmt = sprintf('%%s,X%%0.%sd',num2str(char));
    str = sprintf(fmt,str,i);
end
str = sprintf('%s] = ndgrid(0:dc:c);',str);
eval(str);

X = cell(n,1);
[X{:}] = deal(zeros(size(X1)));
for i = 1:n % transfer X's into a cell array
    eval(sprintf('X{%d} = X%d;',i,i))
end

Pairs = nchoosek(1:n,2);
alpha = cell(n,n);
[alpha{:}] = deal(zeros(size(X1)));
for i = 1:nchoosek(n,2) % create scaled differences
    alpha{Pairs(i,1),Pairs(i,2)} = (X{Pairs(i,1)} .^ 2 - X{Pairs(i,2)} .^ 2) / (sqrt(2) * c);
    alpha{Pairs(i,2),Pairs(i,1)} = -alpha{Pairs(i,1),Pairs(i,2)};
end

x = cell(n,1);
y = cell(n,1);
for i = 1:n % transform to fixed basis
    x{i} = -1 / (2 * n) * (-sqrt(2) * c + 2 * sum(cat(n + 1,alpha{i,:}),n + 1))...
        + 1 / (2 * n) * sqrt(4 * n * X{i} .^ 2 + (-sqrt(2) * c + 2 * sum(cat(n + 1,alpha{i,:}),n + 1)) .^ 2 ...
        - 4 * n * (c ^ 2 / 2 + sum(cat(n + 1,alpha{i,:}) .^ 2,n + 1))); % positive root

    y{i} = -1 / (2 * n) * (-sqrt(2) * c + 2 * sum(cat(n + 1,alpha{i,:}),n + 1))...
        - 1 / (2 * n) * sqrt(4 * n * X{i} .^ 2 + (-sqrt(2) * c + 2 * sum(cat(n + 1,alpha{i,:}),n + 1)) .^ 2 ...
        - 4 * n * (c ^ 2 / 2 + sum(cat(n + 1,alpha{i,:}) .^ 2,n + 1))); % negative root
end

V = cell(n,1);
[V{:}] = deal(zeros(size(X1)));
for i = 1:n % define ith velocity component
    for j = 1:n
        if i < j
            V{i} = V{i} + (-1) ^ (i + j - 1) * X{j};
        elseif i > j
            V{i} = V{i} + (-1) ^ (i + j) * X{j};
        end
    end
end

g = cell(n);
h = cell(n);
for i = 1:n
    for j = 1:n
        g{i,j} = x{j} - States(i,j); % jth component of vector connecting point to ith reference state
        h{i,j} = y{j} - States(i,j); % jth component of vector connecting point to ith reference state
    end
    magnitude1 = sqrt(sum(cat(n + 1,g{i,:}) .^ 2,n + 1));
    magnitude2 = sqrt(sum(cat(n + 1,h{i,:}) .^ 2,n + 1));
    for j = 1:n
        g{i,j} = g{i,j} ./ magnitude1; % convert to components of unit vectors
        h{i,j} = h{i,j} ./ magnitude2; % convert to components of unit vectors
    end
end

v = cell(n,1);
[v{:}] = deal(zeros(size(X1)));
u = cell(n,1);
[u{:}] = deal(zeros(size(X1)));
for i = 1:n
    for j = 1:n
        v{i} = v{i} + V{j} .* g{j,i}; % ith fixed frame velocity due to jth relative frame velocity
        u{i} = u{i} + V{j} .* h{j,i};
    end
    % fix imaginary positions and velocities
    v_temp = v{i};
    x_temp = x{i};
    X_temp = X{i};
    v_temp(imag(x_temp) ~= 0) = 0;
    X_temp(imag(x_temp) ~= 0) = NaN;
    x_temp(imag(x_temp) ~= 0) = 0;
    v{i} = v_temp;
    x{i} = x_temp;
    X{i} = X_temp;

    u_temp = u{i};
    y_temp = y{i};
    u_temp(imag(y_temp) ~= 0) = 0;
    y_temp(imag(y_temp) ~= 0) = 0;
    u{i} = u_temp;
    y{i} = y_temp;
end

H = sum(cat(n + 1,X{:}) .^ 2,n + 1) / 2; %energy function

Phi = arrayfun(@(I) X{Pairs(I,1)} .* X{Pairs(I,2)},1:nchoosek(n,2),'UniformOutput',false);
Phi = sum(cat(n + 1,Phi{:}),n + 1); % potential function

Psi = arrayfun(@(I) (X{Pairs(I,1)} .^ 2 - X{Pairs(I,2)} .^ 2) / 2,1:nchoosek(n,2),'UniformOutput',false);
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
    fprintf('%s, ',c0(1:end - 1))
    fprintf('%s]\n,',c0(end))
end

E_new = arrayfun(@(I) (c0 - States(I,:)) / sqrt(sum((c0 - States(I,:)) .^ 2)),1:n,'UniformOutput',false); % unit vectors toward reference states
X_new = arrayfun(@(I) sqrt(sum((c0 - States(I,:)) .^ 2)),1:n,'UniformOutput',false); % relative frame coordinates
phi_new = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(n,2),'UniformOutput',false);
phi_new = sum(cat(1,phi_new{:}),1); % potential function
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

gamma = zeros(n,npoints); % initialize matrix of sample points
gamma(:,1) = c0'; % first point
phi = zeros(1,npoints);
phi(1) = phi_new;
for k = 2:npoints
    c_new = gamma(:,k - 1)'; % last point
    V_old = V_new;
    for i = 1:n
        c_new = c_new + V_old{i} * step_size; % update position
    end
    gamma(:,k) = c_new';
    E_new = arrayfun(@(I) (c_new - States(I,:)) / sqrt(sum((c_new - States(I,:)) .^ 2)),1:n,'UniformOutput',false); % unit vectors toward reference states
    X_new = arrayfun(@(I) sqrt(sum((c_new - States(I,:)) .^ 2)),1:n,'UniformOutput',false); % relative frame coordinates
    phi_new = arrayfun(@(I) X_new{Pairs(I,1)} .* X_new{Pairs(I,2)},1:nchoosek(n,2),'UniformOutput',false);
    phi_new = sum(cat(1,phi_new{:}),1); % potential function
    phi(k) = phi_new;
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

if n == 2
    if ishandle(1)
        set(0,'CurrentFigure',1)
        cla
    else
        figure(1)
    end
    % Plot trajectory and vector field on total energy surface
    ax = surf(x{1},x{2},H,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.3);
    view(45,15)
    hold on
    surf(y{1},y{2},H,'EdgeColor',get(ax,'EdgeColor'),'FaceColor',get(ax,'FaceColor'),'FaceAlpha',get(ax,'FaceAlpha'))
    plot3(gamma(1,:),gamma(2,:),c ^ 2 / 2 * ones(size(gamma,2)),'LineWidth',2,'Color','k')
    ax = quiver3(x{1},x{2},H,v{1},v{2},zeros(size(H)),'Color','k','LineWidth',1);
    axis([0,c,0,c,0,Inf])
    daspect([1,1,1])
    hold on
    quiver3(y{1},y{2},H,u{1},u{2},zeros(size(H)),'Color',get(ax,'Color'),'LineWidth',get(ax,'LineWidth'))
    
    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.ZLabel.String = '$$H$$';
    ax.ZLabel.Rotation = 0;
    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
    set(ax,'defaulttextinterpreter','latex')
    hold off

    if ishandle(2)
        set(0,'CurrentFigure',2)
        cla
    else
        figure(2)
    end
    % plot trajectory and vector field on potential energy surface
    ax = surf(x{1},x{2},Phi,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.3);
    view(15,45)
    hold on
    surf(y{1},y{2},Phi,'EdgeColor',get(ax,'EdgeColor'),'FaceColor',get(ax,'FaceColor'),'FaceAlpha',get(ax,'FaceAlpha'))
    plot3(gamma(1,:),gamma(2,:),phi,'LineWidth',2,'Color','k')
    ax = quiver3(x{1},x{2},Phi,v{1},v{2},-2 * Psi,'Color','k','LineWidth',1);
    axis([0,c,0,c,0,Inf])
    daspect([1,1,1])
    hold on
    quiver3(y{1},y{2},Phi,u{1},u{2},-2 * Psi,'Color',get(ax,'Color'),'LineWidth',get(ax,'LineWidth'))
    
    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.ZLabel.String = '$$\Phi$$';
    set(ax,'defaulttextinterpreter','latex')
    ax.ZLabel.Rotation = 0;
    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
elseif n == 3
    if ishandle(1)
        set(0,'CurrentFigure',1)
        cla
    else
        figure(1)
    end
    [x1,x2,x3] = ndgrid(0:c/100:c);
    H = (2 * x1 .^ 2 + 2 * x2 .^ 2 + 2 * x3 .^ 2 +...
        (x1 - c / sqrt(2)) .^ 2 + (x2 - c / sqrt(2)) .^ 2 + (x3 - c / sqrt(2)) .^ 2) / 2; % energy function in terms of fixed basis
    tol = 0.01;
    x1 = x1(abs(H - (n - 1) * c ^ 2 / 2) < tol); % find x, y, z coordinates corresponding to set value of H = c ^ 2
    x2 = x2(abs(H - (n - 1) * c ^ 2 / 2) < tol);
    x3 = x3(abs(H - (n - 1) * c ^ 2 / 2) < tol);
    scatter3(x1,x2,x3,'filled','b','o','MarkerEdgeColor','none','MarkerFaceAlpha',0.05)
    view(115,10)
    hold on

    plot3(gamma(1,:),gamma(2,:),gamma(3,:),'LineWidth',3,'Color','k') % level surface of energy function
    hold on
    ax = quiver3(x{1},x{2},x{3},v{1},v{2},v{3},'Color','r','LineWidth',1); % 3D vector field
    hold on
    quiver3(y{1},y{2},y{3},u{1},u{2},u{3},'Color',get(ax,'Color'));
    axis([0,c,0,c,0,c])
    daspect([1,1,1])

    ax = gca;
    ax.XLabel.String = '$$x$$';
    ax.YLabel.String = '$$y$$';
    ax.ZLabel.String = '$$z$$';
    zp = get(get(gca,'ZLabel'),'Position');
    zp(2) = 2 * zp(2); % double z-label position from axis
    set(get(gca,'ZLabel'),'Position',zp)
    ax.ZLabel.Rotation = 0;
    set(ax,'FontSize',16)
    set(ax,'LineWidth',2)
    set(ax,'defaulttextinterpreter','latex')
    hold off

    if save
        exportgraphics(gcf,sprintf('%s/Trajectory n = %d, start = (%0.6f,%0.6f,%0.6d).png',pwd,n,c0(1),c0(2),c0(3)),'Resolution',300)
    end
end




