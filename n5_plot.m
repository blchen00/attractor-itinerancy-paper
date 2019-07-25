% n=5 with & without depression on the same plot
% date 04-24-19
clear

load n5_dep_w40.mat
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
        'bo','MarkerSize',12,'linewidth',1)
    errorbar(stdW,meanNf,err_meanNf,err_meanNf,err_stdW,err_stdW,...
        'bs','MarkerSize',12,'linewidth',1)

    set(gca,'fontsize',fs,'fontname','times','LineWidth',aw)
    xlabel('$\sigma$','FontSize',fs,'Interpreter','latex')
    ylabel('$N_{\rm{states}}$','FontSize',fs,'Interpreter','latex')
    title(strcat('$\mu =',num2str(mu_vec(i_mu)),'$'),'Interpreter','latex')
    box on
    axis square   
    axis([0 0.4 1 32])
%     axis([0 2 1 32])


%     saveas(gcf,strcat('n5_dep9_mean',num2str(mu_vec(i_mu)),'.eps'),'epsc')
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
        'bo','MarkerSize',12,'linewidth',1)

    errorbar(meanW,meanNf,err_meanNf,err_meanNf,err_meanW,err_meanW,...
        'bs','MarkerSize',12,'linewidth',1)

    set(gca,'fontsize',fs,'fontname','times','LineWidth',aw)
    xlabel('$\mu$','FontSize',fs,'Interpreter','latex')
    ylabel('$N_{\rm{states}}$','FontSize',fs,'Interpreter','latex')
    title(strcat('$\sigma =',num2str(sig_vec(i_sig)),'$'),'Interpreter','latex')
    box on
    axis square
    axis([-0.15 0.15 1 32])
%     axis([-0.5 0.5 1 32])


%     saveas(gcf,strcat('n5_dep9_sig',num2str(sig_vec(i_sig)),'.eps'),'epsc')
end

load n5_nodep_w20.mat
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
        'ro','MarkerSize',12,'linewidth',1)
    errorbar(stdW,meanNf,err_meanNf,err_meanNf,err_stdW,err_stdW,...
        'rs','MarkerSize',12,'linewidth',1)

    set(gca,'fontsize',fs,'fontname','times','LineWidth',aw)
    xlabel('$\sigma$','FontSize',fs,'Interpreter','latex')
    ylabel('$N_{\rm{states}}$','FontSize',fs,'Interpreter','latex')
    title(strcat('$\mu =',num2str(mu_vec(i_mu)),'$'),'Interpreter','latex')
    box on
    axis square   
    axis([0 0.4 1 32])
%     axis([0 2 1 32])


%     saveas(gcf,strcat('n5_dep9_mean',num2str(mu_vec(i_mu)),'.eps'),'epsc')
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
        'ro','MarkerSize',12,'linewidth',1)

    errorbar(meanW,meanNf,err_meanNf,err_meanNf,err_meanW,err_meanW,...
        'rs','MarkerSize',12,'linewidth',1)

    set(gca,'fontsize',fs,'fontname','times','LineWidth',aw)
    xlabel('$\mu$','FontSize',fs,'Interpreter','latex')
    ylabel('$N_{\rm{states}}$','FontSize',fs,'Interpreter','latex')
    title(strcat('$\sigma =',num2str(sig_vec(i_sig)),'$'),'Interpreter','latex')
    box on
    axis square
    axis([-0.15 0.15 1 32])
%     axis([-0.5 0.5 1 32])


%     saveas(gcf,strcat('n5_dep9_sig',num2str(sig_vec(i_sig)),'.eps'),'epsc')
end


load n5_nodep_w40.mat
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
        'ko','MarkerSize',12,'linewidth',1)
    errorbar(stdW,meanNf,err_meanNf,err_meanNf,err_stdW,err_stdW,...
        'ks','MarkerSize',12,'linewidth',1)

    set(gca,'fontsize',fs,'fontname','times','LineWidth',aw)
    xlabel('$\sigma$','FontSize',fs,'Interpreter','latex')
    ylabel('$N_{\rm{states}}$','FontSize',fs,'Interpreter','latex')
    title(strcat('$\mu =',num2str(mu_vec(i_mu)),'$'),'Interpreter','latex')
    box on
    axis square   
    axis([0 0.4 1 32])
%     axis([0 2 1 32])


    saveas(gcf,strcat('n5_mean',num2str(mu_vec(i_mu)),'.eps'),'epsc')
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
        'k-o','MarkerSize',12,'linewidth',1)

    errorbar(meanW,meanNf,err_meanNf,err_meanNf,err_meanW,err_meanW,...
        'k-s','MarkerSize',12,'linewidth',1)

    set(gca,'fontsize',fs,'fontname','times','LineWidth',aw)
    xlabel('$\mu$','FontSize',fs,'Interpreter','latex')
    ylabel('$N_{\rm{states}}$','FontSize',fs,'Interpreter','latex')
    title(strcat('$\sigma =',num2str(sig_vec(i_sig)),'$'),'Interpreter','latex')
    box on
    axis square
    axis([-0.15 0.15 1 32])
%     axis([-0.5 0.5 1 32])


    saveas(gcf,strcat('n5_sig',num2str(sig_vec(i_sig)),'.eps'),'epsc')
end