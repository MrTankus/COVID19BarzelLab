%% SEIR where seperating lock-down to areas

clear

global N alpha gamma


%% parameters
N = 200; % num of areas
beta = 1e-3;
alpha = 4e-3;
gamma = 83e-4;
V = 3e3;
K = 1e-2;
pop = 9e6; % Israel population
days = 3e2;
frac_imm = 0.1; % fraction of immigrants 
eta = 1.7e-2;

%%

locs = locs_of_lattice(N); % coordinates of areas
N = size(locs,1);
M = pop/N* max( normrnd(1,0.2,[1,N]),0.1 );

M_imm = M*frac_imm; % number of immigrants from each area

A = build_diff_mat(locs,M_imm);

%% initial conditions
S = 500*rand(1,N);
E = rand(1,N);
I = rand(1,N);
R = rand(1,N);

psi0 = [S;E;I;R];
psi0 = psi0./sum(psi0).*M;

%% Eqs.
% no (I replaced N by M)

der = @(psi,beta) [ -beta.*(psi(1,:).*psi(2,:))./N;
    -gamma*psi(2,:) + beta.*(psi(1,:).*psi(2,:))./N;
    -alpha*psi(3,:) + gamma*psi(2,:);
    alpha*psi(3,:)];

%%

S_days = zeros(days,N);
E_days = zeros(days,N);
I_days = zeros(days,N);
R_days = zeros(days,N);
M_days = zeros(days,N);

psi = psi0;
for id = 1:days
    
    %% color of area
    RG_vec = ( (psi(2,:)+psi(3,:))./M ) < eta;
    
    %% day
    
    beta_vec = beta*RG_vec;
    [psi,M] = go_to_work(A,psi,M,RG_vec);
    if any(psi(:)<0)
        
    end
    t = linspace(8,17,1e3);
    dt = t(2)-t(1);
    for it=1:length(t)
        psi = psi + dt*der(psi,beta_vec);
    end
%     psi = reshape(psi',[4*N,1]);
%     [t,psi] = ode45(@(t,psi) Eq(psi,beta_vec),[8,17],psi);
%     psi = psi(end,:)';
%     psi = reshape(psi,[N,4])';
    

    %% night
    
    [psi,M] = go_home(A,psi,M,RG_vec);
    beta_vec(:) = 0;
    
     if any(psi(:)<0)
        
    end
    t= linspace(17,8+24,1e3);
    dt = t(2)-t(1);
    for it=1:length(t)
        psi = psi + dt*der(psi,beta_vec);
    end
%     psi = reshape(psi',[4*N,1]);
%     [t,psi] = ode45(@(t,psi) Eq(psi,beta_vec),[17,8+24],psi);
%     psi = psi(end,:)';
%     psi = reshape(psi,[N,4])';
    
     if any(psi(:)<0)
        
    end
    
    S_days(id,:) = psi(1,:);
    E_days(id,:) = psi(2,:);
    I_days(id,:) = psi(3,:);
    R_days(id,:) = psi(4,:);
    M_days(id,:) = M;
end

%%

Imax = V/K;

figure
plot(I_days)
xlabel('days')
ylabel('$I_n$','Interpreter','latex')
set(gca,'FontSize',20,'box','on')

figure; hold on
plot(sum(I_days,2),'LineWidth',1)
h = plot(xlim,xlim*0+Imax,'LineWidth',1);
legend(h,{'$I_{max}$'},'Interpreter','latex')
title (['\eta =',num2str(eta) ])
xlabel('days')
ylabel('$I$','Interpreter','latex')
set(gca,'FontSize',20,'box','on')

figure
plot(E_days)
title (['\eta =',num2str(eta) ])
xlabel('days')
ylabel('$E_n$','Interpreter','latex')
set(gca,'FontSize',20,'box','on')

figure
plot(sum(E_days,2))
title (['\eta =',num2str(eta) ])
xlabel('days')
ylabel('$E$','Interpreter','latex')
set(gca,'FontSize',20,'box','on')


figure
plot(M_days)
title (['\eta =',num2str(eta) ])
xlabel('days')
ylabel('$M$','Interpreter','latex')
set(gca,'FontSize',20,'box','on')

figure
plot(sum(M_days,2))
title (['\eta =',num2str(eta) ])
xlabel('days')
ylabel('$N$','Interpreter','latex')
set(gca,'FontSize',20,'box','on')

%%

function der = Eq(psi,beta)
global N alpha gamma
    der = [ -beta'.*(psi(1:end/4).*psi(end/4+1:end/2))./N;
    -gamma*psi(end/4+1:end/2) + beta'.*(psi(1:end/4).*psi(end/4+1:end/2))./N;
    -alpha*psi(end/2+1:3*end/4) + gamma*psi(end/4+1:end/2);
    alpha*psi(end/2+1:3*end/4)];
end
