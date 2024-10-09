%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ERC Project %%%%%%%%%%%%%%%%%%%%%%%
%%%% generating transition probability matrix %%%%%%%%%
%%%% for a repetition (incl. deletion) channel %%%%%%%%
%%%%%%%%% with a max output length limit %%%%%%%%%%%%%%
%%%%%%%%%%%(Poisson, Uniform, Binomial) %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Mohammad Kazemi %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% November 2022 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% inputs: L: channel input length / Rmax: per-bit maximum repetition
%%% output: p: transition probability P(y/x)

%%% the receiver knows maximum of R
%%% reducing out of memory problem
%%% optimized the run-time
%%% limiting repetition per bit to Rmax

function [p,p_Rstar]=TM_RepCh_allR_diffLambda(size_p,I,Rep_pattern,lambda)

RepProb=poisspdf(single(Rep_pattern),lambda);
RepProb=prod(RepProb,2);
sum_prob=sum(RepProb);
RepProb=RepProb/sum_prob; % probability of each repetition pattern
p_Rstar=1-sum_prob;

Np=size(Rep_pattern,1);
Ny_prime=size_p(1); % size of input alphabet
Nx=size_p(2); % size of input alphabet
p=zeros(Ny_prime,Nx);
for j=1:Nx
    for k=1:Np
        ind=(j-1)*Np+k;
%         display(ind);
        p(I(ind),j)=p(I(ind),j)+RepProb(k);
    end
end

p(sum(p,2)==0,:)=[];
