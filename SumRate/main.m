clc;
clear;
close all;
%warning off;
tau_p=2;
Nu = 1; % UE antennas
tc=40;
max_avg = 10; % No of montecarlo
CHA_SETUP = 1;% No of channel setups
T=200;
SNR1=30;
snr1=10.^(SNR1/10);
 
sbl_iter=500;
em_thresh = 0.001;

P=1;
D=128;
N=D;
snr_db0 = -20:10:30;
snr = 10.^(.1*snr_db0);
K = 10;
M = floor(N/3);
% M=N;
L = floor(M/2);

lmmse = 1;
sbl_corr = 1;
sbl_corr_cofem = 1;

SUM_RATE_LMMSE_SIM_LB=zeros(length(snr_db0),1);
SUM_RATE_CORR_SBL_SIM_LB=zeros(length(snr_db0),1);
SUM_RATE_CORR_SBL_COFEM_SIM_LB=zeros(length(snr_db0),1);

LMMSE_NMSE=zeros(length(snr_db0),1);
CORR_SBL_NMSE=zeros(length(snr_db0),1);
CORR_SBL_COFEM_NMSE=zeros(length(snr_db0),1);

for ue_iters = 1:CHA_SETUP
    disp(['CHA_SETUP = ' num2str(ue_iters)]);

    sum_rate_lmmse_sim=zeros(length(snr_db0),1);
    sum_rate_corr_sbl_sim=zeros(length(snr_db0),1);
    sum_rate_corr_sbl_cofem_sim=zeros(length(snr_db0),1);
    
    H=(1/sqrt(N))*exp(1i*2*pi*rand(M,N));
    U=zeros(D,D,K);
    R=0.5;
    for k = 1:K
        %U1=exponentialcorrelation(N,R);  % exponential correlation
        U1 = 0.5*eye(D) + 0.5*ones(D);    % uniform correlation
        U(:,:,k) = U1;
    end
    cov_theta=zeros(D,D,K);

    ind_g=randperm(D,L);
        gamma=zeros(D,1);
        for i=1:L
            gamma(ind_g(i))=1;
        end
        Gamma=diag(gamma);
    for k=1:K
        cov_theta(:,:,k)=(sqrt(Gamma)*U(:,:,k)*sqrt(Gamma));
    end
    A1=dftmtx(D)/sqrt(D);
    x=randperm(D,N);
    A=A1(x,:);
    phi=H*A;

    %% Terms initialization of LMMSE %%
    CCZ_term1=zeros(K,length(snr_db0));
    CCZ_term2=zeros(K,length(snr_db0));
    CCZ_term3=zeros(K,length(snr_db0));

    %% Terms initialization of CORR SBL %%
    CCX_term1=zeros(K,length(snr_db0));
    CCX_term2=zeros(K,length(snr_db0));
    CCX_term3=zeros(K,length(snr_db0));
    
    %% Terms initialization of CORR SBL COFEM %%
    CCY_term1=zeros(K,length(snr_db0));
    CCY_term2=zeros(K,length(snr_db0));
    CCY_term3=zeros(K,length(snr_db0));

    nmse_sbl_corr=zeros(K,max_avg,length(snr_db0));
    nmse_sbl_corr_cofem=zeros(K,max_avg,length(snr_db0));
    nmse_lmmse=zeros(K,max_avg,length(snr_db0));
    nmse_sbl_corr_1=zeros(length(snr_db0),1);
    nmse_sbl_corr_cofem_1=zeros(length(snr_db0),1);
    nmse_lmmse_1=zeros(length(snr_db0),1);

    for mc_iter = 1:max_avg
         fprintf('mc_iter=%d\n',mc_iter);
                h=zeros(N,T*K);
                R_R=zeros(N,N,K);
                R_R1=zeros(N,N,K);
                R_R2=zeros(N,N,K);
                R_avg1 = zeros(N,N);
                R1 = zeros(N,N);
                R_avg = zeros(N,N);
                R_avg2 = zeros(N,N);
               
                for k=1:K
                    [h1,Rh] = channel_generation(D,T,cov_theta(:,:,k),A,N);
                    h(:,1+T*(k-1):T*k)=h1;
                    %Rh=A*cov_theta(:,:,k)*A';
                    
                    R_R(:,:,k)=Rh;
                    R_avg1 = R_avg1+ Rh;
                    R1=R1 + Rh;
                    
                end

                for s = 1:length(snr_db0)  
                     h_hat_corr_sbl1= zeros(N,T*K);
                     h_hat_corr_sbl2= zeros(N,T*K);
                     h_hat_lmmse1 = zeros(N,T*K);
                     y = zeros(M,T*K);
                     R_1 = zeros(N,N);
                     R_2 = zeros(N,N);
                     for k=1:K 
                         noise = sqrt((0.5)/snr1)*(randn(M,T)+1i*randn(M,T));
                         y(:,1+T*(k-1):T*k) = H*h(:,1+T*(k-1):T*k) + noise;
                     end
                     for k=1:K
                         if lmmse == 1
                             [h_hat_lmmse] = LMMSE_ch(R_R(:,:,k),1/snr1,H,y(:,1+T*(k-1):T*k),M);
                             h_hat_lmmse1(:,1+T*(k-1):T*k)=h_hat_lmmse;
                             for t=1:T
                                 nmse_lmmse(k,mc_iter,s) = nmse_lmmse(k,mc_iter,s)+ norm(h_hat_lmmse1(:,1+t*(k-1):t*k)-h(:,1+t*(k-1):t*k))^2/norm(h(:,1+t*(k-1):t*k))^2/T;
                             end
                         end
                     end
                     for k=1:K
                         if sbl_corr == 1
                             [g_hat_corr_sbl,R_g] = SBL_CORR4(U(:,:,k),T,phi,y(:,1+T*(k-1):T*k),1/snr1,sbl_iter,D,M,H);
                             h_hat_corr_sbl = A*g_hat_corr_sbl;
                             R_h=A*R_g*A';
                             R_R1(:,:,k)=R_h;
                             R_avg = R_avg+ R_h;
                             R_1 = R_1 + R_h;
                             h_hat_corr_sbl1(:,1+T*(k-1):T*k)=h_hat_corr_sbl;
                             for t=1:T
                                 nmse_sbl_corr(k,mc_iter,s) = nmse_sbl_corr(k,mc_iter,s)+ norm(h_hat_corr_sbl1(:,1+t*(k-1):t*k)-h(:,1+t*(k-1):t*k))^2/norm(h(:,1+t*(k-1):t*k))^2/T;
                             end
                         end
                     end
                     for k=1:K
                         if sbl_corr_cofem == 1
                             %[g_hat_corr_sbl_cofem,R_g] = SBL_CORR4(U(:,:,k),T,phi,y(:,1+T*(k-1):T*k),1/snr1,sbl_iter,D,M,H);
                             [g_hat_corr_sbl_cofem,R_g,omega_inv_initial,x,omega_g,omega_inv_sum] = SBL_CORR_COFEM(U(:,:,k),T,phi,y(:,1+T*(k-1):T*k),1/snr1,sbl_iter,em_thresh,D,M,H);
                             h_hat_corr_sbl_cofem = A*g_hat_corr_sbl_cofem;
                             R_h=A*R_g*A';
                             R_R2(:,:,k)=R_h;
                             R_avg2 = R_avg2+ R_h;
                             R_2 = R_2 + R_h;
                             h_hat_corr_sbl2(:,1+T*(k-1):T*k)=h_hat_corr_sbl_cofem;
                             for t=1:T
                                 nmse_sbl_corr_cofem(k,mc_iter,s) = nmse_sbl_corr_cofem(k,mc_iter,s)+ norm(h_hat_corr_sbl2(:,1+t*(k-1):t*k)-h(:,1+t*(k-1):t*k))^2/norm(h(:,1+t*(k-1):t*k))^2/T;
                             end
                         end
                     end
                         [Wr_3,Fr_3,Fr_hat_3] = combiner_generation(R_avg,K,N,h_hat_corr_sbl1,M,T);
                         [Wr_4,Fr_4,Fr_hat_4] = combiner_generation(R_avg1,K,N,h_hat_lmmse1,M,T);
                         [Wr_5,Fr_5,Fr_hat_5] = combiner_generation(R_avg,K,N,h_hat_corr_sbl2,M,T);
                    %% SIMULATED LB LMMSE %%
                    for k=1:K
                        CCZ_term1(k,s)=CCZ_term1(k,s) + Wr_4(k,:)*Fr_4*h(:,(k-1)*T+1);
                        for q=1:K
                             CCZ_term2(k,s) = CCZ_term2(k,s) + (abs(Wr_4(k,:)*Fr_4*h(:,(q-1)*T+1)))^2;
                        end
                        CCZ_term3(k,s) = CCZ_term3(k,s) + Wr_4(k,:)*(Fr_4*Fr_4')*Wr_4(k,:)';
                    end
                    %% SIMULATED LB CORR SBL %%
                    for k=1:K
                        CCX_term1(k,s)=CCX_term1(k,s) + Wr_3(k,:)*Fr_3*h(:,(k-1)*T+1);
                        for q=1:K
                             CCX_term2(k,s) = CCX_term2(k,s) + (abs(Wr_3(k,:)*Fr_3*h(:,(q-1)*T+1)))^2;
                        end
                        CCX_term3(k,s) = CCX_term3(k,s) + Wr_3(k,:)*(Fr_3*Fr_3')*Wr_3(k,:)';
                    end
                    %% SIMULATED LB CORR SBL COFEM %%
                    for k=1:K
                        CCY_term1(k,s)=CCY_term1(k,s) + Wr_5(k,:)*Fr_5*h(:,(k-1)*T+1);
                        for q=1:K
                             CCY_term2(k,s) = CCY_term2(k,s) + (abs(Wr_5(k,:)*Fr_5*h(:,(q-1)*T+1)))^2;
                        end
                        CCY_term3(k,s) = CCY_term3(k,s) + Wr_5(k,:)*(Fr_5*Fr_5')*Wr_5(k,:)';
                    end
                 end
    end


   nmse_sbl_corr_2=mean(nmse_sbl_corr,1);
   nmse_lmmse_2=mean(nmse_lmmse,1);
   nmse_sbl_corr_cofem_2=mean(nmse_sbl_corr_cofem,1);
    for s=1:length(snr_db0)
        for mc_iter=1:max_avg
            nmse_sbl_corr_1(s)=nmse_sbl_corr_1(s) + nmse_sbl_corr_2(1,mc_iter,s);
            nmse_sbl_corr_cofem_1(s)=nmse_sbl_corr_cofem_1(s) + nmse_sbl_corr_cofem_2(1,mc_iter,s);
            nmse_lmmse_1(s)=nmse_lmmse_1(s) + nmse_lmmse_2(1,mc_iter,s);
        end
        nmse_sbl_corr_1(s)=nmse_sbl_corr_1(s)/max_avg;
        nmse_sbl_corr_cofem_1(s)=nmse_sbl_corr_cofem_1(s)/max_avg;
        nmse_lmmse_1(s)=nmse_lmmse_1(s)/max_avg;
    end
    CCZ_term1 = CCZ_term1/max_avg;
    CCZ_term2 = CCZ_term2/max_avg;
    CCZ_term3 = CCZ_term3/max_avg;

    CCX_term1 = CCX_term1/max_avg;
    CCX_term2 = CCX_term2/max_avg;
    CCX_term3 = CCX_term3/max_avg;
    
    CCY_term1 = CCY_term1/max_avg;
    CCY_term2 = CCY_term2/max_avg;
    CCY_term3 = CCY_term3/max_avg;    

    %% SIMULATED SE calculation for LMMSE %%
    SINRzz = zeros(K,length(snr_db0));
    SEz1 = zeros(K,length(snr_db0));
    for s = 1:length(snr_db0)
       
        for k=1:K
            SINRzz(k,s)=(abs(CCZ_term1(k,s)))^2/(CCZ_term2(k,s)-((abs(CCZ_term1(k,s)))^2)+(1/snr(s))*CCZ_term3(k,s));
            SEz1(k,s)=(1-(K/tc))*real(log2(1+SINRzz(k,s)));
        end
    end
    %% SIMULATED SE calculation for CORR SBL%%
    SINRxx = zeros(K,length(snr_db0));
    SEx1 = zeros(K,length(snr_db0));
    for s = 1:length(snr_db0)
       
        for k=1:K
            SINRxx(k,s)=(abs(CCX_term1(k,s)))^2/(CCX_term2(k,s)-((abs(CCX_term1(k,s)))^2)+(1/snr(s))*CCX_term3(k,s));
            SEx1(k,s)=(1-(K/tc))*real(log2(1+SINRxx(k,s)));
        end
    end
    %% SIMULATED SE calculation for CORR SBL COFEM%%
    SINRyy = zeros(K,length(snr_db0));
    SEy1 = zeros(K,length(snr_db0));
    for s = 1:length(snr_db0)
       
        for k=1:K
            SINRyy(k,s)=(abs(CCY_term1(k,s)))^2/(CCY_term2(k,s)-((abs(CCY_term1(k,s)))^2)+(1/snr(s))*CCY_term3(k,s));
            SEy1(k,s)=(1-(K/tc))*real(log2(1+SINRyy(k,s)));
        end
    end
   
    %% SUMRATE calculation for LMMSE %%
    for s = 1:length(snr_db0)
        for k=1:K
            sum_rate_lmmse_sim(s)=sum_rate_lmmse_sim(s) + SEz1(k,s);
        end
    end

    %% SUMRATE calculation for CORR SBL%%
    for s = 1:length(snr_db0)
        for k=1:K
            sum_rate_corr_sbl_sim(s)=sum_rate_corr_sbl_sim(s) + SEx1(k,s);
        end
    end
    %% SUMRATE calculation for CORR SBL COFEM%%
    for s = 1:length(snr_db0)
        for k=1:K
            sum_rate_corr_sbl_cofem_sim(s)=sum_rate_corr_sbl_cofem_sim(s) + SEy1(k,s);
        end
    end
    %% 
    SUM_RATE_LMMSE_SIM_LB=SUM_RATE_LMMSE_SIM_LB+sum_rate_lmmse_sim;
    SUM_RATE_CORR_SBL_SIM_LB=SUM_RATE_CORR_SBL_SIM_LB+sum_rate_corr_sbl_sim;
    SUM_RATE_CORR_SBL_COFEM_SIM_LB=SUM_RATE_CORR_SBL_COFEM_SIM_LB+sum_rate_corr_sbl_cofem_sim;

    LMMSE_NMSE=LMMSE_NMSE+nmse_lmmse_1;
    CORR_SBL_NMSE=CORR_SBL_NMSE+nmse_sbl_corr_1;
    CORR_SBL_COFEM_NMSE=CORR_SBL_COFEM_NMSE+nmse_sbl_corr_cofem_1;
end

LMMSE_NMSE=LMMSE_NMSE/CHA_SETUP;
CORR_SBL_NMSE=CORR_SBL_NMSE/CHA_SETUP;
CORR_SBL_COFEM_NMSE=CORR_SBL_COFEM_NMSE/CHA_SETUP;

SUM_RATE_LMMSE_SIM_LB=SUM_RATE_LMMSE_SIM_LB/CHA_SETUP;
SUM_RATE_CORR_SBL_SIM_LB=SUM_RATE_CORR_SBL_SIM_LB/CHA_SETUP;
SUM_RATE_CORR_SBL_COFEM_SIM_LB=SUM_RATE_CORR_SBL_COFEM_SIM_LB/CHA_SETUP;

figure(1)
if lmmse == 1
    semilogy(snr_db0,LMMSE_NMSE,'g-o','linewidth',2,'DisplayName','LMMSE NMSE');
    hold on
end
if sbl_corr == 1
    semilogy(snr_db0,CORR_SBL_NMSE,'c-o','linewidth',2,'DisplayName','SBL CORR NMSE');
    hold on
end
if sbl_corr_cofem == 1
    semilogy(snr_db0,CORR_SBL_COFEM_NMSE,'b-o','linewidth',2,'DisplayName','SBL CORR COFEM NMSE');
    hold on
end
hold off
xlabel('SNR (dB)');
ylabel('NMSE');
legend show
grid on
figure(2)
if lmmse == 1
    plot(snr_db0,SUM_RATE_LMMSE_SIM_LB,'-^r','DisplayName','SIMULATED LB LMMSE');
    hold on
end
if sbl_corr == 1
    plot(snr_db0,SUM_RATE_CORR_SBL_SIM_LB,'--r','DisplayName','SIMULATED LB CORR SBL');
    hold on
end
if sbl_corr_cofem == 1
    plot(snr_db0,SUM_RATE_CORR_SBL_COFEM_SIM_LB,'b--o','DisplayName','SIMULATED LB CORR SBL COFEM');
    hold on
end
hold off
xlabel('SNR (dB)');
ylabel('Sum Rate');
legend show
grid on












































