clear;
close all;
clc;
format long;


% variables that will remain constant in simulations
sbl_corr_cofem_iter = 400;
sbl_corr_iter=400;
sbl_corr_learn_iter=400;
em_thresh=0.001;

% channel setup (Monte carlo loops)
max_avg = 100;

% Select what you want to vary with NMSE 
SNR = 1;
num_L = 0; %number of multipaths
num_T = 0; %number of snapshots
num_D = 0;
num_D_time = 0;
num_M = 0;

%---------
% NMSE vs SNR
if SNR == 1 
    snr_db0 = 15:5:30;
    loop = length(snr_db0);
    N0=128; %128
    D0=128; %128
    M0=70; %60
    K0=10;  %10
    L0=25;  %30
    T0=50; %50
end
% NMSE vs L
if num_L == 1
    L0 = 20:5:55;
    loop = length(L0);
    snr_db0 = 25;
    N0=128;
    D0=128;
    M0=60;
    K0=10;
    T0=100;
end
% NMSE vs T
if num_T == 1
     T0 =[1,5,10,50,100,200,400,600,800,1000];
%     T0 =[1,10,100];
    loop = length(T0);
    snr_db0 = 10;
    N0=128;
    D0=128;
    M0=64;
    L0=32;
    K0=8;
end
% NMSE vs D
if num_D == 1
    D0 = 80:20:220;
    %D0 = [256,512,1024];
    loop = length(D0);
    snr_db0 = 25;
    N0 = D0;
    M0 = 60; %60
    L0 = 30; %30
    K0 = 10;
    T0 = 50;
end
% Time vs D
if num_D_time == 1
    D0 = 2.^[8:12];
    loop = length(D0);
    snr_db0 = 25;
    N0 = D0;
    M0 = D0/2; %60
    L0 = round(M0/3); %30 D/4
    K0 = 10;
    T0 = 50;
end
% NMSE vs M
if num_M ==1
    M0 = 48:10:128;
    loop = length(M0);
    snr_db0 = 25;
    D0 = 128;
    N0 = D0;
    L0 = 30; %30
    K0 = 10;
    T0 = 50;
    %30 - 128
end
% initialization of error metrics

nmse_sbl_1 = zeros(loop,1);
nmse_sbl_corr_1 = zeros(loop,1);
nmse_lmmse_1=zeros(loop,1);
NMSE_SBL=zeros(loop,1);
NMSE_SBL_CORR=zeros(loop,1);
NMSE_LMMSE=zeros(loop,1);
TIME_SBL_CORR = zeros(loop,1);
TIME_SBL_CORR_COFEM = zeros(loop,1);
TIME_LMMSE = zeros(loop,1);
%-------------------------------------
%=====================================
% select the algorithms you want to run
sbl_corr_cofem = 1;
sbl_corr = 1;
lmmse=1;
%=======================================
%---------------------------------------

wt=waitbar(0,'Initializing...');
curr_iter=0;
for loop_iter = 1:loop
    nmse_1=0;
    nmse_2=0;
    nmse_4=0;
    time_1=0;
    time_2=0;
    time_4=0;
    nmse_sblcorr_cofem=zeros(max_avg,1);
    nmse_sbl_corr=zeros(max_avg,1);
    nmse_lmmse=zeros(max_avg,1);
    time_sbl_corr_cofem = zeros(max_avg,1);
    time_sbl_corr = zeros(max_avg,1);
    time_lmmse = zeros(max_avg,1);

    for mc_iter = 1:max_avg
%         fprintf('mc_iter=%d\n',mc_iter);
        if SNR == 1
            snr = 10.^(snr_db0(loop_iter)/10);
            N=N0;
            D=D0;
            M=M0;
            K=K0;
            L=L0;
            T=T0;
        end
        if num_L == 1
            L=L0(loop_iter);
            snr = 10.^(snr_db0/10);
            N=N0;
            D=D0;
            M=M0;
            K=K0;
            T=T0;
        end
        if num_T == 1
            T=T0(loop_iter);
            snr = 10.^(snr_db0/10);
            N=N0;
            D=D0;
            M=M0;
            K=K0;
            L=L0;
        end
        if num_D == 1
            D = D0(loop_iter);
            snr = 10.^(snr_db0/10);
            N=N0(loop_iter);
            M=M0;
            K=K0;
            L=L0;
            T = T0;
        end
        if num_D_time == 1
            D = D0(loop_iter);
            snr = 10.^(snr_db0/10);
            N=N0(loop_iter);
            M=M0(loop_iter);
            K=K0;
            L=L0(loop_iter);
            T = T0;
        end
        if num_M == 1
            M = M0(loop_iter);
            snr = 10.^(snr_db0/10);
            N=N0;
            D=D0;
            K=K0;
            L=L0;
            T = T0;
        end
        U = 0.5*eye(D) + 0.5*ones(D);
        [h,G, ind_g, array_response,cov_theta] = channel_generation(D,N,U,L,T);
               
        [phi,W] = dictionary_generation(N,M, array_response);
        
        g_hat_sbl=zeros(D,T);
        g_hat_sblcorr=zeros(D,T);
        h_hat_sbl_cofem=zeros(D,T);
        h_hat_sblcorr=zeros(D,T);
        h_hat_lmmse = zeros(D,T);

        
        %data generation
        noise = sqrt(0.5/snr)*(randn(M,T)+1i*randn(M,T));
        y = phi*G.' + noise;
       
        if sbl_corr_cofem == 1
             U = 0.5*eye(D) + 0.5*ones(D);
             c1 = inv(U).*(cov_theta.');
             k=10;

             tic;
             [g_hat_sblcorr_cofem,R_g1,x2,x3,omega_g1,omega_inv_sum,sig1,omega_c1] = SBL_CORR_2(U,T,phi,y,1/snr,sbl_corr_cofem_iter,em_thresh,D,M);
             time_sbl_corr_cofem(mc_iter) = toc;
                
             h_hat_sbl_cofem = array_response*g_hat_sblcorr_cofem;
             for i=1:T
                 nmse_sblcorr_cofem(mc_iter) = nmse_sblcorr_cofem(mc_iter)+ norm(h_hat_sbl_cofem(:,i)-h(:,i))^2/norm(h(:,i),2)^2/T;
             end

             nmse_1=nmse_1 + nmse_sblcorr_cofem(mc_iter);
             time_1 = time_1 + time_sbl_corr_cofem(mc_iter);
             %if isnan(nmse_sblcorr_cofem(mc_iter))
               %  break;
             %end
        end
        
        if sbl_corr == 1
            U = 0.5*eye(D) + 0.5*ones(D);
            tic;
            [g_hat_sblcorr,R_g,x,omega_g,omega_c,sigm] = SBL_CORR_1(U,T,phi,y,1/snr,sbl_corr_iter,em_thresh,D,M,L,cov_theta);
            time_sbl_corr(mc_iter) = toc;

            h_hat_sblcorr = array_response*g_hat_sblcorr;
            for i=1:T
            nmse_sbl_corr(mc_iter) = nmse_sbl_corr(mc_iter)+ norm(h_hat_sblcorr(:,i)-h(:,i))^2/norm(h(:,i))^2/T;
            end
            nmse_2=nmse_2 + nmse_sbl_corr(mc_iter);
            time_2 = time_2 + time_sbl_corr(mc_iter);
        end
        if lmmse == 1
            
            tic;
            Rh=array_response*cov_theta*array_response';
            [h_hat_lmmse] = LMMSE(Rh,1/snr,W,y,M);
            time_lmmse(mc_iter) = toc;

            for i=1:T
            nmse_lmmse(mc_iter) = nmse_lmmse(mc_iter)+ norm(h_hat_lmmse(:,i)-h(:,i))^2/norm(h(:,i))^2/T;
            end
            nmse_4=nmse_4 + nmse_lmmse(mc_iter);
            time_4 = time_4 + time_lmmse(mc_iter);
        end
        curr_iter = 1+curr_iter;
        waitbar(curr_iter/(max_avg*loop),wt,sprintf('%0.1f%% done',curr_iter/(max_avg*loop)*100))
    end


    nmse_sbl_1(loop_iter) = nmse_1/(max_avg);
    NMSE_SBL(loop_iter)=20*log10(nmse_sbl_1(loop_iter));
    nmse_sbl_corr_1(loop_iter) = nmse_2/(max_avg);
    NMSE_SBL_CORR(loop_iter)=20*log10(nmse_sbl_corr_1(loop_iter));
    nmse_lmmse_1(loop_iter) = nmse_4/(max_avg);
    NMSE_LMMSE(loop_iter) = 20*log10(nmse_lmmse_1(loop_iter));
    TIME_SBL_CORR_COFEM(loop_iter) = time_1/(max_avg);
    TIME_SBL_CORR(loop_iter) = time_2/(max_avg);
    TIME_LMMSE(loop_iter) = time_4/(max_avg);
    
    
    if SNR == 1
        if sbl_corr_cofem == 1
            semilogy(snr_db0,nmse_sbl_1,'g-o','linewidth',2,'DisplayName','SBL CORR COFEM');
            hold on
        end
        if sbl_corr == 1
            semilogy(snr_db0,nmse_sbl_corr_1,'y-x','linewidth',2,'DisplayName','SBL CORR');
            hold on
        end
        if lmmse == 1
            semilogy(snr_db0,nmse_lmmse_1,'r-*','linewidth',1,'DisplayName','LMMSE');
            hold on
        end
        hold off
        xlabel('SNR (dB)');
        ylabel('NMSE');
    end
    
    if num_L == 1
        if sbl_corr_cofem == 1
            semilogy(L0,nmse_sbl_1,'g-o','linewidth',2,'DisplayName','SBL CORR COFEM');
            hold on
        end
        if sbl_corr == 1
            plot(L0,nmse_sbl_corr_1,'y-x','linewidth',2,'DisplayName','SBL CORR');
            hold on
        end
        if lmmse == 1
            semilogy(L0,nmse_lmmse_1,'r-*','linewidth',1,'DisplayName','LMMSE');
            hold on
        end
        hold off
        %xlabel('Number of Multipaths L');
        xlabel('sparsity');
        ylabel('NMSE');
        %ylim([10^(-3) inf]);
    end
     
     if num_T == 1
         if sbl_corr_cofem == 1
            semilogx(T0,NMSE_SBL,'b-o','linewidth',2,'DisplayName','SBL');
            hold on
        end
        if sbl_corr == 1
            semilogx(T0,NMSE_SBL_CORR,'y-x','linewidth',2,'DisplayName','SBL CORR');
            hold on
        end
        if lmmse == 1
            semilogx(T0,NMSE_LMMSE,'r-*','linewidth',1,'DisplayName','LMMSE');
            hold on
        end
        hold off
        xlabel('Number of Snapshots T');
        ylabel('NMSE');
     end

     if num_D == 1
        if sbl_corr_cofem == 1
            semilogy(D0,nmse_sbl_1,'g-o','linewidth',2,'DisplayName','SBL CORR COFEM');
            hold on
        end
        if sbl_corr == 1
            semilogy(D0,nmse_sbl_corr_1,'y-x','linewidth',2,'DisplayName','SBL CORR');
            hold on
        end
        if lmmse == 1
            semilogy(D0,nmse_lmmse_1,'r-*','linewidth',1,'DisplayName','LMMSE');
            hold on
        end
        hold off
        xlabel('D');
        ylabel('NMSE');
     end
     
     if num_D_time == 1
        if sbl_corr_cofem == 1
            semilogx(D0,TIME_SBL_CORR_COFEM,'g-o','linewidth',2,'DisplayName','SBL CORR COFEM');
            hold on
        end
        if sbl_corr == 1
            semilogx(D0,TIME_SBL_CORR,'y-x','linewidth',2,'DisplayName','SBL CORR');
            hold on
        end
        if lmmse == 0
            semilogx(D0,TIME_LMMSE,'r-*','linewidth',1,'DisplayName','LMMSE');
            hold on
        end
        grid on;
        hold off
        xlabel('D');
        ylabel('Time Taken');
        legend('SBL COFEM','SBL CORR','LMMSE');
        xticklab = cellstr(num2str(round(log2(D0(:))),'2^{%d}'));
        set(gca,'XTick',D0,'XTickLabel',xticklab,'TickLabelInterpreter','tex');
     end

     if num_M == 1
        if sbl_corr_cofem == 1
            semilogy(M0,nmse_sbl_1,'g-o','linewidth',2,'DisplayName','SBL CORR COFEM');
            hold on
        end
        if sbl_corr == 1
            plot(M0,nmse_sbl_corr_1,'y-x','linewidth',2,'DisplayName','SBL CORR');
            hold on
        end
        if lmmse == 1
            semilogy(M0,nmse_lmmse_1,'r-*','linewidth',1,'DisplayName','LMMSE');
            hold on
        end
        hold off
        %xlabel('Number of Multipaths L');
        xlabel('No of obervations (M)');
        ylabel('NMSE');
    end
    
    legend show
    grid on
end

 delete(wt)



















