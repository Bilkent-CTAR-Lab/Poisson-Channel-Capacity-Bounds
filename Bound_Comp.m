%%% Poisson-repeat channel LB/UB comparison
%%% LB/UB-proposed (the receiver just knows Rmaxout)

L=5;
Rmax=3;
% R=floor(L/3);
R=5;
% R=Rmax*L/2;
tic;

% eps=1e-5; d=[eps,.01:.01:.19,.2:.1:.9,1-eps];
% eps=1e-5; d=.85,.9,.95,1-eps;
% step=.01;d=.9:step:1

d=linspace(0,1,11);
% d=[0, .9 ,1];
lambda=-log(d);

% lambda=.1:.1:.4
% d=exp(-lambda)

[p_tot,I_tot,Rep_pattern_tot,R_out_tot,p_Rstar_tot]=Transition_Matrix_RepCh_allRv3_par(L,Rmax,R,1);
% p_tot=p_tot;I_tot=I_tot;Rep_pattern_tot=Rep_pattern_tot;R_out_tot=R_out_tot;p_Rstar_tot=p_Rstar_tot;
% p_tot=p_tot;I_tot=I_tot;Rep_pattern_tot=Rep_pattern_tot;R_out_tot=R_out_tot;

size_p=size(p_tot);
toc

tic
LB=zeros(1,length(d));
UB=zeros(1,length(d));
capp=zeros(1,length(d));
H=zeros(1,length(d));
p_Rstar=zeros(1,length(d));
for j=1:length(d)
    %     j
    if d(j)==0
        LB(j)=1;
        UB(j)=1;
    elseif d(j)==1
        LB(j)=0;
        UB(j)=0;
    else
        [p,p_Rstar(j)]=TM_RepCh_allR_diffLambda(size_p,I_tot,Rep_pattern_tot,lambda(j));
        %     p_Rstar=1-poisscdf(R,L*lambda(j));
        %     p_Rstar=(1-poisscdf(Rmax,lambda(j))^L)+1-poisscdf(R,L*lambda(j)); % union probability
        capp(j)=BAA(p,.05);
        co=1e5;poiss_pdf=poisspdf(0:max(co*lambda(j),co),L*lambda(j));
        H(j)=-sum(poiss_pdf(poiss_pdf>0).*log2(poiss_pdf(poiss_pdf>0)));
        %         H2=-p_Rstar(j)*log2(p_Rstar(j))-(1-p_Rstar(j))*log2((1-p_Rstar(j)));

        LB(j)=max(0,(1-p_Rstar(j))*capp(j)-H(j))/L;
        UB(j)=(1-p_Rstar(j))*capp(j)/L+p_Rstar(j);
    end
end
UB_Cheraghchi=[1,.61,.5,.41,.335,.275,.205,.16,.095,.045,0]; % UB-Cheraghchi 0:.1:1
LB_Drinea=[1,.73,.57,.44,.35,.29,.23,.18,.15,.12,.1,.085,.075,.065,.05,.04,.025,.017,.01,.005,0];
%     step=.05;d:0:.05:1
UB_Cheraghchi2=interp1(0:.1:1,UB_Cheraghchi,d); % UB-Cheraghchi 0:step:1

% plot(d,LB,d,UB,'LineWidth',2)
% plot(d,LB,d,UB,d,lambda/9,'-.',[d(1:end-1),1],UB_Cheraghchi2(end-length(d)+1:end),'--','LineWidth',2)
% plot(lambda,LB,lambda,UB,lambda,lambda/9,'-.',lambda,UB_Cheraghchi2,'--','LineWidth',2)
%     plot(d,LB./lambda,d,UB./lambda,d,ones(1,length(lambda))/9,'-.',[d(1:end-1),1],UB_Cheraghchi2(end-length(d)+1:end)./lambda,'--','LineWidth',2)

% plot(d,LB,d,UB,d,lambda/9,'-.',d,UB_Cheraghchi2,'--','LineWidth',2)
% legend('LB-proposed','UB-proposed','LB-1/9','UB-Cheraghchi')

plot(d,UB,'LineWidth',2)

% xlabel('$\frac{C_{\lambda}}{\lambda}$')
% xlabel('d (e^{-\lambda})')
xlabel('e^{-\lambda}')
% ylabel('bounds - bits per channel use')
ylabel('bits per channel use')
% ylabel('normalized bounds (over $\lambda$)')
grid on; set(gcf,'color','w'); fontsize(gca,20,'points')
hold on
%     L
%     lambda
LB./lambda
p_Rstar
% capp
% H
% (1-p_Rstar).*capp-H
(UB_Cheraghchi2-UB)./UB_Cheraghchi2

toc

% save data_9_3_9.mat -v7.3;