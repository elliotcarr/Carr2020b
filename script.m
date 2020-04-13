clc, close all, clear all

%% Setup

% save_figs = true;
save_figs = false;

% Plotting uses the export_fig package (download from here: https://github.com/altmany/export_fig)
font_size = 34;
line_width = 3;
if save_figs
    path_name = '../../Paper/Figures/'; % folder to save figures
    addpath('../export_fig-master') % path to export_fig
end

% Figure 1
Case = 'A'; delta = 5e-2; plot_ts = true; table_ts = false; kindx = [15]; tspan = [1e-2,1];
% Case = 'B'; delta = 5e-2; plot_ts = true; table_ts = false; kindx = [15]; tspan = [1e-2,1];

% Table 1
% Case = 'A'; delta = 1e-2; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = [];
% Case = 'A'; delta = 1e-4; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = [];
% Case = 'A'; delta = 1e-6; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = [];
% Case = 'B'; delta = 1e-2; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = [];
% Case = 'B'; delta = 1e-4; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = [];
% Case = 'B'; delta = 1e-6; plot_ts = false; table_ts = true; kindx = [1,5,10,15]; tspan = [];

N = 501; % number of nodes
m = 15; % maximum number of moments
L = 1.0; % length of medium
cb = 1.0; % concentration of species 1 at x = 0
n = 3; % number of species
D = @(x) 0.1 + 0.05*sin(10*x); % diffusivity
Lbnd = {cb,0,0}; % left boundary data.
Rbnd = {0,0,0}; % right boundary data.
cint = zeros(N,n); % initial condition
if strcmp(Case,'A')
    mu = [0.8,0.4,0.1]; % reaction rates
elseif strcmp(Case,'B')
    mu = [0.8,0.4,0.35]; % reaction rates 
end

%% Steady-state solution
f1 = zeros(N,1);
f2 = zeros(N,2);
[x,A,b,Mt,map] = discretisation(n,D,mu,L,N,Lbnd,Rbnd);
cinf = A\b;
cinf = reshape(cinf,N,n);
cinf = full(cinf);

%% Moments
M = cell(m+1,1);
Mbar = cell(m+1,1);
LbndM = Lbnd; LbndM{1} = 0; % left boundary data
RbndM = Rbnd; % right boundary data

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
pt = N; % evaluate transition time at node "pt".
tau = zeros(m,2);
for k = 1:m
    for i = 1:n
        tau(k,i) = (M{k+1}(pt,i)/(k*M{k}(pt,i)))*log((M{k+1}(pt,i)...
            /(factorial(k)*delta))*(k*M{k}(pt,i)/M{k+1}(pt,i))^k);
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
    
    colors = [228,26,28; 55,126,184; 77,175,74]/255;
    str = cell(n,1);
    for i = 1:n
        str{i} = num2str(i);
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
        leg = legend(str,'Interpreter','LaTeX','Fontsize',font_size-6,'Location','NorthEast');
        leg.Title.String = 'Species';
        set(gca,'Fontsize',font_size-2,'TickLabelInterpreter','latex','Xtick',xticks,...
            'XMinorTick','on','YTick',[0,0.5,1],'YMinorTick','on','TickDir','out',...
            'Color',0.97*ones(1,3))
        xlabel('$x$','Interpreter','LaTeX','Fontsize',font_size)
        ylabel('$c_{i}(x,t)$','Interpreter','LaTeX','Fontsize',font_size)
        text(0.35*L,0.8,['Case ',Case'],'Interpreter','LaTeX','Fontsize',font_size-2)
        text(0.35*L,0.7,['$t = ',num2str(tspan(p)),'$'],'Interpreter','LaTeX','Fontsize',font_size-2)
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
        plot(x,cn(:,indxp(i),i),'--o','Color',colors(i,:),'LineWidth',line_width,...
            'MarkerIndices',round(linspace(1,N,20)),'MarkerSize',8,'MarkerFaceColor',colors(i,:));
    end
    axis([0 L 0 1])
    leg = legend(ha,str,'Interpreter','LaTeX','Fontsize',font_size-6,'Location','NorthEast');
    leg.Title.String = 'Species';
    set(gca,'Fontsize',font_size-2,'TickLabelInterpreter','latex','Xtick',xticks,'XMinorTick','on',...
        'YTick',[0,0.5,1],'YMinorTick','on','TickDir','out','Color',0.97*ones(1,3))
    xlabel('$x$','Interpreter','LaTeX','Fontsize',font_size)
    ylabel('$c_{i}(x,t)$','Interpreter','LaTeX','Fontsize',font_size)
    text(0.35*L,0.8,['Case ',Case],'Interpreter','LaTeX','Fontsize',font_size-3)
    text(0.35*L,0.7,'$t\rightarrow\infty$','Interpreter','LaTeX','Fontsize',font_size-3)
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