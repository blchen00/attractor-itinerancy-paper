% Random nets of bistable units without synaptic depression
% date: 04-18-19
% test N=5 without depression 
% parameters are wii = 15, theat = 5;
% initial state index = 9
% update: 04-20-19
% improve the searching for attractors
% update:04-23-19
% try to redo n=5 with depression
% update: 04-24-19
% fixed fxt searching; combining non-interacting guess solutions with
% random initial samplings



clear
fs = 32; lw=2.5;  aw=1.5;
sympref('HeavisideAtOrigin', 1);  %  heaviside(0)=1
opts = odeset('RelTol',1e-8,'InitialStep',1e-6);


FLAG_DEP = 1;               % 0 = no depression; 1 = with depression
attractor_digit = 9;        % initial state (01001)
% attractor_digit = 1;        % initial state (00001)
N = 5;

taur = 0.01*ones(N,1);
taus = 0.05*ones(N,1);
taud = 0.25*ones(N,1);
ts = taus./taur;    % dimensionless synaptic time constant
td = taud./taur;    % dimensionless depression time constant

alpha = 1./ts;
beta = 1./td;

rmax = 50*ones(N,1);
p0 = 0.5*ones(N,1);
alp = ones(N,1);

a = FLAG_DEP.*p0.*taud.*rmax;         % dimensionless depression parameter
b = alp.*p0.*taus.*rmax;    % dimensionless synpatic parameter

% location of the cusp point
wc = 4*(a+b+1)./b;
thc = 2+log(a+b+1);

% choose self-excitation and threshold in the bistable region
Wii = 40.*eye(N);
% Wii = 20.*eye(N);       % reduce self-excitation to 20 when depression is off
Theta = 5.*ones(N,1);

% std & mean for Wij
sig_vec = [linspace(0,0.2,9),0.3:0.1:2];
mu_vec = -0.6:0.2:0.6;

sig_vec = [0.001,0.01,0.1,0.2,0.3];  %  parameters in old code using fxpt solver
mu_vec = [-0.1,0,0.1];

% sig_vec = [0.001,0.01];  %  parameters in old code using fxpt solver
% mu_vec = [-0.1,0,0.1];
% mu_vec = 0;
% sig_vec = 0;
% 
% sig_vec = [0.001,0.01,0.1,0.2,0.3,0.4:0.2:2];  %  parameters in old code using fxpt solver
% mu_vec = -0.4:0.1:0.4;

Nsig = length(sig_vec);
Nmu = length(mu_vec);
Nsample = 20;  


Wij_data = nan.*ones(N,N,Nsample,Nsig,Nmu);
stdW_data = nan.*ones(Nsample,Nsig,Nmu);
meanW_data = nan.*ones(Nsample,Nsig,Nmu);

Na_data = nan.*ones(Nsample,Nsig,Nmu);
Nf_data = nan.*ones(Nsample,Nsig,Nmu);

ra_data = nan.*ones(N,2^N,Nsample,Nsig,Nmu);
rf_data = nan.*ones(N,2^N,Nsample,Nsig,Nmu);

ra_bin_data = nan.*ones(N,2^N,Nsample,Nsig,Nmu);
rf_bin_data = nan.*ones(N,2^N,Nsample,Nsig,Nmu);

% load n5_nodepression_init9_part.mat
% tmp_counters = [i_mu,i_sig,i_sample];

% define trial solutions to speed up fxpt searching
r0_trial = zeros(N,2^N);
r0_trial_bin = dec2bin(0:1:2^N-1);
for ii=1:2^N
    for jj=1:N
        r0_trial(jj,ii)=bin2dec(r0_trial_bin(ii,jj));
    end
end

sinf =@(x) b.*x./(1+(a+b).*x);   % steady sol for s(r)
dinf =@(x) 1./(1+a.*x);          % steady sol for d(r)

r0_trial = 0.6.*FLAG_DEP.*r0_trial;
s0_trial = sinf(r0_trial(:,:));
d0_trial = dinf(r0_trial(:,:));
Xinit_vec = [r0_trial; s0_trial; d0_trial];  % trial fxpts


for i_mu=1:Nmu
    mu = mu_vec(i_mu);
    for i_sig=1:Nsig
        sig = sig_vec(i_sig);
        for i_sample=1:Nsample
            Wij = zeros(N);
            Wij_vec = [];
            for ii=1:N-1
                tmp_triu = sig.*randn(N-ii,1)+mu;
                tmp_tril = sig.*randn(N-ii,1)+mu;
                Wij = Wij + diag(tmp_triu,ii) + diag(tmp_tril,-ii);
                Wij_vec = [Wij_vec;tmp_triu;tmp_tril];
            end
%             % make sum of cross-connections zero to each unit
%             Wij = Wij - sum(Wij,2)*ones(1,N)/(N-1) - mu; 
%             Wij = Wij - diag(diag(Wij));
%             for i = 1:N; 
%                 for j = i+1:N; Wij_vec = [Wij_vec;Wij(i,j)]; end
%                 for j =i-1:-1:1; Wij_vec = [Wij_vec;Wij(i,j)]; end
%             end
            
            W = Wii + Wij;
            Wij_data(:,:,i_sample,i_sig,i_mu) = Wij(:,:);
            
            stdWij = std(Wij_vec);
            meanWij = mean(Wij_vec);
            stdW_data(i_sample,i_sig,i_mu) = stdWij;
            meanW_data(i_sample,i_sig,i_mu) = meanWij;


            %% finding stable fixed points 
            ti = 0; t0 = 10; tf = 1000;     
            tau = 0;
            Iapp = 0;

            Xa_vec = [];                % attractors array
            Na = 0;                    % number of attractors
            Ntrials = 0;
            
            tic
            while (Na<2^N && Ntrials<3000 )

                switch FLAG_DEP
                    case 1
                        if(Ntrials<2^N)
                            Xinit = Xinit_vec(:,Ntrials+1);
                        else
                        Xinit = rand(3*N,1);
                        end
                    case 0
                        if(Ntrials<2^N)
                            Xinit = Xinit_vec(:,Ntrials+1);
                        else
                        
                        rinit = rand(N,1);
                        sinit = zeros(N,1);
                        dinit = ones(N,1);   
                        Xinit = [rinit;sinit;dinit];
                        end
                end
                  

                sol0 = ode45(@(t,X) ...
                    odefun(t,X,alpha,beta,a,b,W,Theta,Iapp,t0,tau,FLAG_DEP),...
                    [ti tf], Xinit,opts);
                X0 = sol0.y(:,end);

                new_states = 1;             % Assume it is a new attractor   
                for i = 1:Na                % Test all prior states
                    if (norm(X0 - Xa_vec(:,i)) < 1e-3 ) % if it is close to a prior one
                        new_states = 0;                 % then not a new state
                        Ncounts(i) = Ncounts(i) + 1;    % Record number of instances 
                    end
                end
                if (new_states == 1)        % If it is a new state
                    Na = Na + 1;            % Add 1 to state numbers
                    Xa_vec(:,Na) = X0;      % Record the state vector
                    Ncounts(Na) = 1;        % Add 1 instance of this state
                end

                Ntrials = Ntrials + 1;
            end
            toc
            
            % binarize the state vector
            ra_vec = Xa_vec(1:N,:);
            ra_bin = ra_vec>0.5;
            ra_bin_digit = transpose(bin2dec(num2str(transpose(ra_bin))));
            
            % if initial state is not in found attractors
            % then skip to the next sampled network
            if ~ismember(attractor_digit,ra_bin_digit) 
                disp('initial state is not in found attractors')
                continue;                              
            end                                   
            
            % if some states occur more than once in found attractors
            % then faulty solutions, skip to the next sampled network
            if length(unique(ra_bin_digit)) < length(ra_bin_digit)
                disp('some states occur more than once in found attractors')
                continue;
            end                                       
            
            % estimate the probability of taking a stable state
            Pstates = Ncounts/Ntrials; 
            Estates = -log(Pstates);
            disp([Na, Ntrials])
            disp(ra_vec)
            disp(ra_bin)

            % save Nattractors to Nattractors_data array
            Na_data(i_sample,i_sig,i_mu) = Na;
            for i_states=1:Na
                ra_data(:,i_states,i_sample,i_sig,i_mu) = ...
                    ra_vec(:,i_states);
                
                ra_bin_data(:,i_states,i_sample,i_sig,i_mu) = ...
                    ra_bin(:,i_states);
            end

            %% perturb an attractor state 
            ntrial = 32;
            Iapp_vec = linspace(0,5,ntrial); 
            tau_vec = linspace(1,200,ntrial); 

            % duration and amplitude mesh
            [taum,Iappm] = meshgrid(tau_vec,Iapp_vec);
            zm = zeros(size(taum));
            ngrid = numel(taum);
            
            % start with the same attractor specified as attractor_digit
            Xa_idx = find(ra_bin_digit==attractor_digit);
            Xa = Xa_vec(:,Xa_idx);

            parfor itrial=1:ngrid
                Iapp = Iappm(itrial);             
                tau = taum(itrial);

                sol = ode45(@(t,X) ...
                    odefun(t,X,alpha,beta,a,b,W,Theta,Iapp,t0,tau,FLAG_DEP),...
                    [ti tf], Xa, opts);
                Xf = sol.y(:,end);
                
                [~,Xf_idx] = min(vecnorm(Xf - Xa_vec));
                rf_bin_digit = ra_bin_digit(Xf_idx);
                
%                 rf = Xf(1:N);
%                 rf_bin = rf>0.5;
%                 rf_bin_digit = transpose(bin2dec(num2str(transpose(rf_bin))));

                zm(itrial) = rf_bin_digit;
            end 
            
            rf_digits = unique(zm);
            rf_idxs = find(ismember(ra_bin_digit,rf_digits));
            rf_vec = ra_vec(:,rf_idxs);
            rf_bin = rf_vec>0.5;
            

            Nf = length(rf_digits);
            Nf_data(i_sample,i_sig,i_mu) = Nf;

            for i_final = 1:Nf
                rf_data(:,i_final,i_sample,i_sig,i_mu)=...
                    rf_vec(:,i_final);
                
                rf_bin_data(:,i_final,i_sample,i_sig,i_mu)=...
                    rf_bin(:,i_final);
            end

            disp([i_mu,i_sig,i_sample,stdWij,meanWij,Nf,Na]);
            
        end % end of sampling random networks
        
%         filename = strcat('n',num2str(N),'_mean0_sig',num2str(sig),...
%             '_mu_',num2str(mu),'.mat');
%         save(filename)
        
    end % end of sigma loop
    
    

end % end of mu loop

save('n5_depression_init9.mat');

for i_mu = 1:length(mu_vec)
    figure(i_mu)
    hold on
    meanW = mean(meanW_data(:,:,i_mu));
    std_meanW = std(meanW_data(:,:,i_mu));

    stdW = mean(stdW_data(:,:,i_mu));
    std_stdW = std(stdW_data(:,:,i_mu));

    meanNa = nanmean(Na_data(:,:,i_mu));
    std_meanNa = nanstd(Na_data(:,:,i_mu));

    meanNf = nanmean(Nf_data(:,:,i_mu));
    std_meanNf = nanstd(Nf_data(:,:,i_mu));
    
    err_stdW = std_stdW./2;
    err_meanNa = std_meanNa./2;
    err_meanNf = std_meanNf./2;

    errorbar(stdW,meanNa,err_meanNa,err_meanNa,err_stdW,err_stdW,...
        'o','MarkerSize',12,'linewidth',1)
    errorbar(stdW,meanNf,err_meanNf,err_meanNf,err_stdW,err_stdW,...
        'o','MarkerSize',12,'linewidth',1)

    set(gca,'fontsize',fs,'fontname','times','LineWidth',aw)
    xlabel('$\sigma$','FontSize',fs,'Interpreter','latex')
    ylabel('$N_{\rm{states}}$','FontSize',fs,'Interpreter','latex')
    title(strcat('$\mu =',num2str(mu_vec(i_mu)),'$'),'Interpreter','latex')
    box on
    axis square   
    axis([0 0.4 1 32])
%     axis([0 2 1 32])


    saveas(gcf,strcat('n5_dep9_mean',num2str(mu_vec(i_mu)),'.eps'),'epsc')
end

% permute data sets
meanW_data_new = permute(meanW_data,[1,3,2]);
stdW_data_new = permute(stdW_data,[1,3,2]);
Na_data_new = permute(Na_data,[1,3,2]);
Nf_data_new = permute(Nf_data,[1,3,2]);
ra_data_new = permute(ra_data,[1,2,3,5,4]);
rf_data_new = permute(rf_data,[1,2,3,5,4]);

for i_sig = 1:length(sig_vec)
    figure(i_sig+length(mu_vec))
    hold on
    meanW = mean(meanW_data_new(:,:,i_sig));
    std_meanW = std(meanW_data_new(:,:,i_sig));

    stdW = mean(stdW_data_new(:,:,i_sig));
    std_stdW = std(stdW_data_new(:,:,i_sig));

    meanNa = nanmean(Na_data_new(:,:,i_sig));
    std_meanNa = nanstd(Na_data_new(:,:,i_sig));

    meanNf = nanmean(Nf_data_new(:,:,i_sig));
    std_meanNf = nanstd(Nf_data_new(:,:,i_sig));

    err_meanW = std_meanW./2;
    err_meanNa = std_meanNa./2;
    err_meanNf = std_meanNf./2;
    
    errorbar(meanW,meanNa,err_meanNa,err_meanNa,err_meanW,err_meanW,...
        'o','MarkerSize',12,'linewidth',1)

    errorbar(meanW,meanNf,err_meanNf,err_meanNf,err_meanW,err_meanW,...
        's','MarkerSize',12,'linewidth',1)

    set(gca,'fontsize',fs,'fontname','times','LineWidth',aw)
    xlabel('$\mu$','FontSize',fs,'Interpreter','latex')
    ylabel('$N_{\rm{states}}$','FontSize',fs,'Interpreter','latex')
    title(strcat('$\sigma =',num2str(sig_vec(i_sig)),'$'),'Interpreter','latex')
    box on
    axis square
    axis([-0.15 0.15 1 32])
%     axis([-0.5 0.5 1 32])


    saveas(gcf,strcat('n5_dep9_sig',num2str(sig_vec(i_sig)),'.eps'),'epsc')
end
        



function dX = odefun(t,X,alpha,beta,a,b,Wij,Theta,Iapp,t0,tau,FLAG_DEP)
    n = numel(X)/3;
    r = X(1:n,:);
    s = X((1:n)+n,:);
    d = X((1:n)+2*n,:);
    
    dr = -r + F(Wij*s - Theta + Iext(t,t0,tau,Iapp));
    ds = alpha.*(-s + b.*r.*d.*(1-s));
    dd = FLAG_DEP.*beta.*(1 - d - a.*r.*d);

    dX = [dr;ds;dd];
end

function y = F(x)
    y = 1./(1 + exp(-x));
end
 
function stimuli = Iext(t,t0,tau,Iapp)
    stimuli = Iapp.*(heaviside(t-t0)-heaviside(t-t0-tau));
end



% % plot the phase diagram 
% figure
% basin_plt = pcolor(taum,Iappm,zm);
% % colormap(lines(Nattrators))  
% set(basin_plt,'edgecolor','none')
% caxis([0 Nattrators-1])
% colorbar('Ticks',0:Nattrators-1)
% axis square
% set(gca,'fontsize',fs,'fontname','times')
% xlabel('$\tau_{\rm{dur}}$','FontSize',fs,'Interpreter','latex')
% ylabel('$I_{\rm{app}}$','FontSize',fs,'Interpreter','latex')