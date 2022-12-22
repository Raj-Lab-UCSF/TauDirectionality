%% 0. Loading data
clear; clc; close all; rng(0);
matdir = '/Users/justintorok/Documents/MATLAB/TauDirectionality_ChrisPaper_Project/TauDirectionality/raw_data_mouse';
figdir = '/Users/justintorok/Documents/MATLAB/TauDirectionality_ChrisPaper_Project/Figures';
load([matdir filesep 'eNDM_mousedata.mat'],'Networks','tpts','netlabs','data426','seed426','labnms');
C_ret = Networks.ret; C_ant = Networks.ant; C_nd = Networks.nd;
C_ret = (C_ret - diag(diag(C_ret)));
C_ant = (C_ant - diag(diag(C_ant)));
C_nd = (C_nd - diag(diag(C_nd)));

taudata_all = struct;
tpts_all = struct;
seed_all = struct;
labnms_all = labnms;
for i = 1:length(labnms)
    taudata_all.(labnms{i}) = data426.(labnms{i});
    tpts_all.(labnms{i}) = tpts.(labnms{i});
    seed_all.(labnms{i}) = seed426.(labnms{i});
end
clearvars -except C_ret C_ant C_nd taudata_all tpts_all labnms_all seed_all netlabs matdir figdir   

load([matdir filesep 'KaufmanDiamond_datasets_dat&seed.mat'],'tpts','data426','seed426','labnms');
labnms_all = [labnms_all.' labnms];
for i = 1:length(labnms)
    taudata_all.(labnms{i}) = data426.(labnms{i});
    tpts_all.(labnms{i}) = tpts.(labnms{i});
    seed_all.(labnms{i}) = seed426.(labnms{i});
end
clearvars -except C_ret C_ant C_nd taudata_all tpts_all labnms_all seed_all netlabs matdir figdir

%% 1. Graph analysis
% G_nd = digraph(C_ret);
% deg_Gnd = centrality(G_nd,'outdegree','Importance',G_nd.Edges.Weight);
% C_nd = logical(C_nd);
% L_ret = eye(426) - pinv(diag(sum(C_ret,2)))*C_ret;
% deg = diag(sum(C_ret,1));
% L_ret = eye(426) - (deg^(-1/2) * C_ret * deg^(-1/2));
% L_ret = eye(426) - deg^(-0.5)*C_ret*deg^(-0.5);
% L_ret = L_ret.';
% L_ret2 = genLplcns(C_ret);
% [v,d] = eig(L_ret); d = real(diag(d)); v = real(v);
% [v2,d2] = eig(L_ret2); 
% d2 = abs(diag(d2)); v2 = abs(v2);
% [dsort,sortinds] = sort(d); 
% vsort = v(:,sortinds);
% [d2sort,sortinds2] = sort(d2);
% v2sort = v2(:,sortinds2);
L_ret = genLplcns(C_ret);
[v_ret,d_ret] = eig(L_ret); d_ret = abs(diag(d_ret)); v_ret = abs(v_ret);
[dretsort,sortinds] = sort(d_ret); vretsort = v_ret(:,sortinds);
L_ant = genLplcns(C_ant);
[v_ant,d_ant] = eig(L_ant); d_ant = abs(diag(d_ant)); v_ant = abs(v_ant);
[dantsort,sortinds] = sort(d_ant); vantsort = v_ant(:,sortinds);

%% 2. Figures
hurtado_end = taudata_all.IbaHippInj(:,end);
hurtado_end_nonan = hurtado_end(~isnan(hurtado_end));
outdeg = sum(C_ret,2); outdeg_hurtado = outdeg(~isnan(hurtado_end));
indeg = sum(C_ret,1).'; indeg_hurtado = indeg(~isnan(hurtado_end));
u1_ret = v_ret(:,1); u1_ret_hurtado = u1_ret(~isnan(hurtado_end));
u1_ant = v_ant(:,1); u1_ant_hurtado = u1_ant(~isnan(hurtado_end));

figure; 
subplot(2,2,1);
scatter(outdeg_hurtado,hurtado_end_nonan,'bo','filled'); lsline;
xlabel('Out Degree'); ylabel('Hurtado End')
legend(sprintf('R = %.2f',corr(outdeg_hurtado,hurtado_end_nonan)));
title('Out Degree vs. Hurtado');
set(gca,'FontSize',16);

subplot(2,2,2);
scatter(indeg_hurtado,hurtado_end_nonan,'ro','filled'); lsline;
xlabel('In Degree'); ylabel('Hurtado End')
legend(sprintf('R = %.2f',corr(indeg_hurtado,hurtado_end_nonan)));
title('In Degree vs. Hurtado');
set(gca,'FontSize',16);

subplot(2,2,3);
scatter(u1_ret_hurtado,hurtado_end_nonan,'bd','filled'); lsline;
xlabel('u1_r_e_t'); ylabel('Hurtado End')
legend(sprintf('R = %.2f',corr(u1_ret_hurtado,hurtado_end_nonan)));
title('u1_r_e_t vs. Hurtado');
set(gca,'FontSize',16);

subplot(2,2,4);
scatter(u1_ant_hurtado,hurtado_end_nonan,'rd','filled'); lsline;
xlabel('u1_a_n_t'); ylabel('Hurtado End')
legend(sprintf('R = %.2f',corr(u1_ant_hurtado,hurtado_end_nonan)));
title('u1_a_n_t vs. Hurtado');
set(gca,'FontSize',16);

%% 3. Conn. from seed
ibahippinj_seed = logical(seed_all.IbaHippInj);
ibahippinj_end = taudata_all.IbaHippInj(:,2);
ibahippinj_end(logical(ibahippinj_seed)) = NaN;
ibahippinj_end_nonan = ibahippinj_end(~isnan(ibahippinj_end));
Cout_ibaseed = C_ret(ibahippinj_seed,:).';
Cin_ibaseed = C_ret(:,ibahippinj_seed);
Cout_ibaseed = Cout_ibaseed(~isnan(ibahippinj_end));
Cin_ibaseed = Cin_ibaseed(~isnan(ibahippinj_end));

figure; 
subplot(1,2,1);
scatter(Cout_ibaseed,ibahippinj_end_nonan,'bo','filled'); lsline;
xlabel('Conn. from CA3 Seed'); ylabel('IbaHippInj End')
legend(sprintf('R = %.2f',corr(Cout_ibaseed,ibahippinj_end_nonan)));
title('Conn. from Seed vs. IbaHippInj');
set(gca,'FontSize',16);

subplot(1,2,2);
scatter(Cin_ibaseed,ibahippinj_end_nonan,'ro','filled'); lsline;
xlabel('Conn. to CA3 Seed'); ylabel('IbaHippInj End')
legend(sprintf('R = %.2f',corr(Cin_ibaseed,ibahippinj_end_nonan)));
title('Conn. to Seed vs. IbaHippInj');
set(gca,'FontSize',16);

%% 4. S vs. Aggregation
% load([cd filesep 'SampleFiles' filesep 'beta_gamma_curve_finerange.mat']);
% figure('Units','inches','Position',[0 0 10 9]); hold on;
% cmap = hsv(size(biasmat,2));
% smat = 0.5 - 0.5*biasmat;
% legcell = cell(1,length(fraclist));
% for i = 1:length(fraclist)
%     legcell{i} = sprintf('f = %.1f',fraclist(i));
% end
% for i = 1:size(biasmat,2)
%     gbratio = gammalist(8:24)/beta;
%     plot(gbratio,smat(8:24,i,end),'x-','Color',cmap(i,:),'LineWidth',3);
% end
% xlabel('\gamma/\beta'); ylabel('s'); title('Steady State Bias vs. Aggregation');
% legend(legcell,'Location','southeast');
% set(gca,'FontSize',20,'xscale','log');

%% S1. Functions
    function L = genLplcns(mat)
    
        Dr = sum(mat,2);
        Dc = sum(mat,1);
        small = find(Dr < 0.05 * mean(Dr));
        Dr(small(:)) = 0.05 * mean(Dr);
        small = find(Dc < 0.05 * mean(Dc));
        Dc(small(:)) = 0.05 * mean(Dc);
        Dr = diag(Dr);
        Dc = diag(Dc);
        
        L = eye(size(mat)) - ((Dr^-(1/2)) * mat * (Dc^-(1/2)));
    end

