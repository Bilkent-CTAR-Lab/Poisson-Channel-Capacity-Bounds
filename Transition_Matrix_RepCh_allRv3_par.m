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
%%% output: p: transition probability P(y/x), I: transistion matrix's
%%% pattern, Rep_pattern: all possible repetition patterns, R_out: output
%%% lengths, p_star: probability that the output is not in the considered
%%% outputs of TM

%%% the receiver knows maximum of R
%%% reducing out of memory problem
%%% optimized the run-time
%%% limiting repetition per bit to Rmax
%%% using appropriate classes (int8,...)
%%% added p_star as output

function [p,I,Rep_pattern,R_out,p_Rstar]=Transition_Matrix_RepCh_allRv3_par(L,Rmax,R,lambda)

Nx=2^L; % size of input alphabet

if L<=4
    Rep_pattern=dec2base(0:((Rmax+1)^L)-1,Rmax+1,L)-'0'; % all possible repetition patterns
    Rep_pattern(Rep_pattern>9)=Rep_pattern(Rep_pattern>9)-7; % fixing the issue that MATLAB considers A as 17 not 10
    I=sum(Rep_pattern,2)>R;
    Rep_pattern(I,:)=[]; % removing the rows that leads to output lengths of more than R
else % to resolve the 'out of memory' problem
    sup=4;
    step_size=(Rmax+1)^sup;
    N_step=(Rmax+1)^(L-sup);
    Rep_pattern=zeros(0,'int8');
%     Rep_pattern_cell=cell(1,N_step);
    parfor j=1:N_step
        Rep_pattern_int=dec2base((j-1)*step_size:j*step_size-1,Rmax+1,L)-'0'; % all possible repetition patterns (intermediate)
        Rep_pattern_int=int8(Rep_pattern_int);
        Rep_pattern_int(Rep_pattern_int>9)=Rep_pattern_int(Rep_pattern_int>9)-7; % fixing the issue that MATLAB considers A as 17 not 10
        I=sum(Rep_pattern_int,2)>R;
        Rep_pattern_int(I,:)=[]; % removing the rows that leads to output lengths of more than R
%         int=(j-1)*step_size+1:j*step_size;
%         Rep_pattern_cell{j}=Rep_pattern_int;
    Rep_pattern=[Rep_pattern;Rep_pattern_int];
    end
%     Rep_pattern=[];
%     for j=1:N_step
%         Rep_pattern=[Rep_pattern;Rep_pattern_cell{j}];
%     end
%     Rep_pattern(sum(abs(Rep_pattern),2)==0)=[];
end
Rep_pattern=int8(Rep_pattern);
Np=size(Rep_pattern,1);

% R=L*Rmax;

if lambda<=1e2
    RepProb=poisspdf(single(Rep_pattern),lambda);
    RepProb=prod(RepProb,2);
    sum_prob=sum(RepProb);
    RepProb=RepProb/sum_prob; % probability of each repetition pattern
else % to resolve NAN due to high lambda
    alpha=zeros(1,R+1);
    parfor j=0:R
        alpha(j+1)=exp((j-R)*log(lambda)+sum(log(j+1:R)));
    end
    alpha=alpha/sum(alpha);
    RepProb=alpha(Rep_pattern+1);
    sum_prob=sum(RepProb);
    RepProb=RepProb/sum_prob; % probability of each repetition pattern
end
RepProb=single(RepProb);
p_Rstar=1-sum_prob;

X=2*(dec2bin(0:Nx-1,L)-'0')-1; % input alphabet (binary string)
X=int8(X);

Y=zeros(0,'int8');
parfor x=1:Nx
    Y_nonbinary=Rep_pattern.*repmat(X(x,:),Np,1);
    Y_temp=zeros(Np,R); % polar form of Y_nonbinary
    for j=1:Np
        temp=zeros(0,'int8');
        for k=1:L
            temp=[temp,sign(Y_nonbinary(j,k))*ones(1,abs(Y_nonbinary(j,k)),'int8')];
        end
        Y_temp(j,1:length(temp))=temp;
    end
    Y_temp=unique(Y_temp,'stable','rows');
    Y=[Y;Y_temp];
end
% Y0=Y;
% Ny0=size(Y,1); % output size after pruning

Y=unique(Y,'stable','rows'); % removing duplicate rows (pruning)
R_out=sum(abs(Y),2);
Ny_prime=size(Y,1); % output size after pruning

% [~,I]=ismember(Y0,Y,'rows');

p=zeros(Ny_prime,Nx);
% I_temp=zeros(Np,1);
I=zeros(Np,Nx);
parfor jj=1:Nx
    Y_nonbinary=Rep_pattern.*repmat(X(jj,:),Np,1);
    Y_temp=zeros(Np,R); % polar form of Y_nonbinary
    for j=1:Np
        temp=[];
        for k=1:L
            temp=[temp,sign(Y_nonbinary(j,k))*ones(1,abs(Y_nonbinary(j,k)),'int8')];
        end
        Y_temp(j,1:length(temp))=temp;
    end 
    [~,I_temp]=ismember(Y_temp,Y,'rows');
    I(:,jj)=I_temp;
    p_temp=zeros(Ny_prime,1,'single');
    for k=1:Np
        p_temp(I_temp(k))=p_temp(I_temp(k))+RepProb(k);
    end
    p(:,jj)=p_temp;
end
I=int32(I(:));

% Y=[];
% parfor x=1:Nx
%     Y_nonbinary=Rep_pattern.*repmat(X(x,:),Np,1);
%     Y_temp=zeros(Np,R); % polar form of Y_nonbinary
%     for j=1:Np
%         temp=[];
%         for k=1:L
%             temp=[temp,sign(Y_nonbinary(j,k))*ones(1,abs(Y_nonbinary(j,k)))];
%         end
%         Y_temp(j,1:length(temp))=temp;
%     end
%     Y_temp=unique(Y_temp,'stable','rows');
%     Y=[Y;Y_temp];
% end
% Y0=Y;
% % Ny0=size(Y,1); % output size after pruning
% 
% Y=unique(Y,'stable','rows'); % removing duplicate rows (pruning)
% R_out=sum(abs(Y),2);
% Ny_prime=size(Y,1); % output size after pruning
% 
% [~,I]=ismember(Y0,Y,'rows');
% 
% p=zeros(Ny_prime,Nx);
% parfor j=1:Nx
%     p_temp=zeros(Ny_prime,1);
%     for k=1:Np
%         ind=(j-1)*Np+k;
%         p_temp(I(ind))=p_temp(I(ind))+RepProb(k);
%     end
%     p(:,j)=p_temp;
% end