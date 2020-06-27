clc, close all, clear all

%% Setup

% save_figs = true;
save_figs = false;
colors = [161,218,180; 65,182,196; 44,127,184; 37,52,148]/255;
colors = colors(end:-1:1,:);

% Plotting uses the export_fig package (download from here: https://github.com/altmany/export_fig)
font_size = 34;
line_width = 3;
if save_figs
    path_name = '../../Paper revision/Figures/'; % folder to save figures
    addpath('../export_fig-master') % path to export_fig
end

% Figure 1
Case = 'A'; delta = 1e-2; plot_ts = true; table_ts = false; kindx = [15]; tspan = [1e-2,1]; spatial = false;
% Case = 'B'; delta = 1e-2; plot_ts = true; table_ts = false; kindx = [15]; tspan = [1e-2,1]; spatial = false;

% Figure 2
% Case = 'C'; delta = 1e-2; plot_ts = true; table_ts = false; kindx = [15];tspan = [1e-2,2]; spatial = false;

% Figure 3
% Case = 'C'; delta = 1e-6; plot_ts = true; table_ts = false; kindx = [15]; tspan = [1e-2,1]; spatial = true;

% Table 1
% Case = 'A'; delta = 1e-2; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = []; spatial = false;
% Case = 'A'; delta = 1e-4; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = []; spatial = false;
% Case = 'A'; delta = 1e-6; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = []; spatial = false;
% Case = 'B'; delta = 1e-2; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = []; spatial = false;
% Case = 'B'; delta = 1e-4; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = []; spatial = false;
% Case = 'B'; delta = 1e-6; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = []; spatial = false;

% Table 2
% Case = 'C'; delta = 1e-2; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = []; spatial = false;
% Case = 'C'; delta = 1e-4; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = []; spatial = false;
% Case = 'C'; delta = 1e-6; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = []; spatial = false;

N = 501; % number of nodes
m = 15; % maximum number of moments
L = 1.0; % length of medium
if strcmp(Case,'A') 
    n = 3; % number of species
    D = @(x) 0.1 + 0.05*sin(10*x); % diffusivity
    mu = [-0.8 0 0; 0.8 -0.4 0; 0 0.4 -0.1]; % reaction rates
    cb = 1.0; % concentration of species 1 at x = 0
    cint = zeros(N,n); % initial condition
    Lbnd = {cb,0,0}; % left boundary data.
elseif strcmp(Case,'B')
    n = 3; % number of species
    D = @(x) 0.1 + 0.05*sin(10*x); % diffusivity
    mu = [-0.8 0 0; 0.8 -0.4 0; 0 0.4 -0.35]; % reaction rates 
    cb = 1.0; % concentration of species 1 at x = 0
    cint = zeros(N,n); % initial condition
    Lbnd = {cb,0,0}; % left boundary data.
elseif strcmp(Case,'C')
    n = 4; % number of species
    %mu = [-0.8 0 0; 0.8*0.5 -0.4 0; 0.8*0.5 0.4*0.75 -0.1]; % reaction rates
    D = @(x) 1.0*(x<0.35) + 0.01*(x>0.35 && x < 0.65) + 1.0*(x>0.65); % diffusivity
    mu = repmat([0.4 0.8 0.2 1.2],4,1) .* [-1 0.5 0.25 0.25; 0 -1 0.3 0.7; 0 0 -1 1; 0.1 0.9 0 -1]';
    cb = 1.0; % concentration of species 1 at x = 0
    cint = zeros(N,n); % initial condition
    Lbnd = {1.0,0.5,0.25,0.7}; % left boundary data. 
end

%% Steady-state solution
f1 = zeros(N,1);
f2 = zeros(N,2);
[x,A,b,Mt,map] = discretisation(n,D,mu,L,N,Lbnd);
cinf = A\b;
cinf = reshape(cinf,N,n);
cinf = full(cinf);

%% Moments
M = cell(m+1,1);
Mbar = cell(m+1,1);
LbndM = Lbnd; LbndM{1} = 0; % left boundary data

Mvec = zeros(N*n,1);
% Moment 0
for i = 1:n
    Mvec(map(1:N,i)) = cinf(1:N,i) - cint(1:N,i);
end
for i = 1:n
    Mbar{1}(:,i) = Mvec(map(1:N,i));
end
% Moments 1,...,m
for k = 1:m+1
    f = k*Mvec;
    f(map(1,1:n)) = 0;
    Mvec = A\f;
    for i = 1:n
        Mbar{k+1}(:,i) = Mvec(map(1:N,i));
    end
end
for k = 0:m
    for i = 1:n
        M{k+1}(:,i) = Mbar{k+1}(:,i) ./ (cinf(:,i) - cint(:,i));
    end
end

%% Transition times
if spatial
    
    tau = zeros(N,n);
    for pt = 1:N % evaluate transition time at node "pt".
        k = m;
        for i = 1:n
            tau(pt,i) = (M{k+1}(pt,i)/(k*M{k}(pt,i)))*log((M{k+1}(pt,i)...
                /(factorial(k)*delta))*(k*M{k}(pt,i)/M{k+1}(pt,i))^k);
        end
    end
    
    figure;
    set(gcf,'Color','w')
    str = cell(n,1);
    for i = 1:n
        str{i} = ['Species ',num2str(i)];
    end
    xticks = 0:0.5:1;
    tau(1,:) = zeros(1,n);
    for i = 1:n
        plot(x,tau(:,i),'Color',colors(i,:),'LineWidth',2)
        hold on
    end
    leg = legend(str,'Interpreter','LaTeX','Fontsize',font_size-16,'Location','SouthEast');
    set(gca,'Fontsize',font_size-12,'TickLabelInterpreter','latex','Xtick',xticks,...
        'XMinorTick','on','YTick',[0,100,200,250],'YMinorTick','on','TickDir','out',...
        'Color','w')
    xlabel('$x$','Interpreter','LaTeX','Fontsize',font_size-10)
    ylabel('$\tau_{i}(x)$','Interpreter','LaTeX','Fontsize',font_size-10)
    drawnow
    if save_figs
        feval('export_fig',[path_name,'Case',Case,'_Spatial'],'-pdf')
    end
    return
    
else
    
    pt = N; % evaluate transition time at node "pt".
    tau = zeros(m,n);
    for k = 1:m
        for i = 1:n
            tau(k,i) = (M{k+1}(pt,i)/(k*M{k}(pt,i)))*log((M{k+1}(pt,i)...
                /(factorial(k)*delta))*(k*M{k}(pt,i)/M{k+1}(pt,i))^k);
        end
    end
    
end

%% Solve transient problem
options = odeset('Mass',Mt,'MassSingular','yes','RelTol',1e-11,'AbsTol',1e-8,'MaxStep',1e-2,...
    'MStateDependence','none','JPattern',A,'BDF','on');

tau = tau(kindx,:);
[tspan,indx] = sort([reshape(tau',1,length(kindx)*n),tspan]);
indxp = zeros(length(kindx)*n,1);
for k = 1:length(kindx)
    for i = 1:n
        indxp((k-1)*n+i) = find(indx == ((k-1)*n+i));
    end
end

F = @(t,c) -A*c+b;
[~,cnt] = ode15s(F,[0,tspan],reshape(cint,n*N,1),options);
cnt = cnt(2:end,:); % remove initial condition
cnt = cnt'; 
   
cn = zeros(N,length(tspan),2);
for i = 1:n
    if length(tspan) == 1
        cntp = cnt(:,end);
        cn(:,1,i) = cntp((i-1)*N+1:i*N);
    else
        for p = 1:length(tspan)
            cntp = cnt(:,p);
            cn(:,p,i) = cntp((i-1)*N+1:i*N);
        end
    end
end

%% Plot transient solution at transition time estimates
if plot_ts
    
    str = cell(n,1);
    for i = 1:n
        str{i} = ['Species ',num2str(i)];
    end    
    xticks = 0:0.5:1;
    
    % Plot species concentrations at early times specified in original tspan vector.
    for p = 1:(length(tspan)-n)
        figure;
        vec = get(gcf,'Position');
        set(gcf,'Color','w','Position',vec.*[1,1,1,1])
        if ~ismember(indx(p),[1:n])
            for i = 1:n
                plot(x,cn(:,p,i),'-','LineWidth',line_width,'Color',colors(i,:))
                hold on
            end
        end
        axis([0 L 0 1])
        if p == 1
            leg = legend(str,'Interpreter','LaTeX','Fontsize',font_size-6,'Location','NorthEast');
        end
        %leg.Title.String = 'Species';
        set(gca,'Fontsize',font_size-2,'TickLabelInterpreter','latex','Xtick',xticks,...
            'XMinorTick','on','YTick',[0,0.5,1],'YMinorTick','on','TickDir','out',...
            'Color','w')
        xlabel('$x$','Interpreter','LaTeX','Fontsize',font_size)
        ylabel('$c_{i}(x,t)$','Interpreter','LaTeX','Fontsize',font_size)
        title(['Case ',Case,' [$t = ',num2str(tspan(p)),'$]'],'Interpreter','LaTeX','Fontsize',font_size-3)
        if strcmp(Case,'A') || strcmp(Case,'C')
            if p == 1
                text(-0.15*L,-0.25,'(a)','Interpreter','LaTeX','Fontsize',font_size)
            else
                text(-0.15*L,-0.25,'(b)','Interpreter','LaTeX','Fontsize',font_size)
            end
        else
            if p == 1
                text(-0.15*L,-0.25,'(d)','Interpreter','LaTeX','Fontsize',font_size)
            else
                text(-0.15*L,-0.25,'(e)','Interpreter','LaTeX','Fontsize',font_size)
            end
        end
        drawnow
        if save_figs
            feval('export_fig',[path_name,'Case',Case,'_Time',num2str(p)],'-pdf')
        end
    end
    
    % Plot each species concentration at its finite transition time.
    figure;
    vec = get(gcf,'Position');
    set(gcf,'Color','w','Position',vec.*[1,1,1.0,1])
    hold on
    ha = zeros(n,1);
    for i = 1:n
        ha(i) = plot(x,cinf(:,i),'-','Color',colors(i,:),'LineWidth',line_width);
    end
    for i = 1:n
        plot(x,cn(:,indxp(i),i),'-o','Color',colors(i,:),'LineWidth',1,...
            'MarkerIndices',round(linspace(1,N,20)),'MarkerSize',12,'MarkerFaceColor',colors(i,:));
    end
    axis([0 L 0 1])
    set(gca,'Fontsize',font_size-2,'TickLabelInterpreter','latex','Xtick',xticks,'XMinorTick','on',...
        'YTick',[0,0.5,1],'YMinorTick','on','TickDir','out','Color','w')
    xlabel('$x$','Interpreter','LaTeX','Fontsize',font_size)
    ylabel('$c_{i}(x,t)$','Interpreter','LaTeX','Fontsize',font_size)
    title(['Case ',Case,' [$t\rightarrow\infty$]'],'Interpreter','LaTeX','Fontsize',font_size-3)
    if strcmp(Case,'A') || strcmp(Case,'C')
        text(-0.15*L,-0.25,'(c)','Interpreter','LaTeX','Fontsize',font_size-1)
    else
        text(-0.15*L,-0.25,'(f)','Interpreter','LaTeX','Fontsize',font_size-1)
    end
    box on
    drawnow
    
    if save_figs
        feval('export_fig',[path_name,'Case',Case,'_steady_state'],'-pdf')
    end
    
end

%% Prints rows of Table 1
if table_ts
    for k = 1:length(kindx)
        for i = 1:n
            ik = (k-1)*n + i;
            if i == 1
                if k == 1
                    fprintf('$10^{-%i}$ & %i & %1.2f & \\num{%1.4e} & ',-log10(delta),kindx(k),...
                        tau(k,i),(cn(pt,indxp(ik),i)-cinf(pt,i))/(cint(pt,i)-cinf(pt,i)))
                else
                    fprintf('& %i & %1.2f & \\num{%1.4e} & ',kindx(k),tau(k,i),...
                        (cn(pt,indxp(ik),i)-cinf(pt,i))/(cint(pt,i)-cinf(pt,i)))
                end
            elseif i > 1 && i < n
                fprintf('%1.2f & \\num{%1.4e} & ',tau(k,i),...
                    (cn(pt,indxp(ik),i)-cinf(pt,i))/(cint(pt,i)-cinf(pt,i)))
            else
                fprintf('%1.2f & \\num{%1.4e}\\\\\n',tau(k,i),...
                    (cn(pt,indxp(ik),i)-cinf(pt,i))/(cint(pt,i)-cinf(pt,i)))
            end
        end
    end
end