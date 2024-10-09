%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% ERC Project %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Blahut-Arimoto Algorithm (BAA) %%%%%%%%%%%%%
%%%%%%%% numerically computes DMC capacity %%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Mohammad Kazemi %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% October 2022 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% input: p: transition probability P(y/x)
%%% output: C: capacity

function [C,Q]=BAA_par(p,epsilon)

if nargin==1
    epsilon=0.005; % stopping criterion (capacity precision)
end
Nx=size(p,2); % size of input alphabet
Ny=size(p,1); % size of output alphabet

%% initialization

% prodd=0;
% while prodd==0
%     Q = rand(1,Nx); 
%     prodd=prod(Q);
% end
% Q = Q/sum(Q); % initial input distribution
Q = ones(1,Nx)/Nx; % initial input distribution

Phi = zeros(Nx,Ny); % Phi(x,y) 

%% BAA

Delta_C=Nx;
while Delta_C>epsilon

    parfor y=1:Ny % updating Phi
        p_temp=p(y,:);
        p_sum=0;
        for x_prime=1:Nx
            p_sum = p_sum+Q(x_prime)*p_temp(x_prime);
        end
        for x=1:Nx
            Phi(x,y) = Q(x)*p_temp(x)/p_sum;
        end
    end
    
    r=ones(1,Nx);
    parfor x=1:Nx % updating r
        p_temp=p(:,x);
        for y=1:Ny
            r(x) = r(x)*(Phi(x,y)^p_temp(y));
        end
    end
    C=log2(sum(r));
    c=r./Q;
    Delta_C=max(log2(c))-C;
    Q = r/sum(r);
end

% Delta_C_min=0.0005; % stoppinf criterion
% Nx=size(p,2); % size of input alphabet
% Ny=size(p,1); % size of output alphabet
% 
% %% initialization
% 
% r = rand(1,Nx); 
% r = r/sum(r); % initial input distribution
% 
% q = zeros(Nx,Ny); % q(x,y) 
% 
% %% BAA
% 
% Delta_C=1;
% C=0;
% while Delta_C>Delta_C_min
% 
%     for y=1:Ny % updating q
%         p_sum=0;
%         for x_prime=1:Nx
%             p_sum = p_sum+r(x_prime)*p(y,x_prime);
%         end
%         for x=1:Nx
%             q(x,y) = r(x)*p(y,x)/p_sum;
%         end
%     end
%     
%     r=ones(1,Nx);
%     for x=1:Nx % updating r
%         for y=1:Ny
%             r(x) = r(x)*(q(x,y)^p(y,x));
%         end
%     end
%     r = r/sum(r);
%     
%     C_old=C;
%     C=0;
%     for x=1:Nx
%         for y=1:Ny
%             if q(x,y)~=0
%                 C = C + r(x)*p(y,x)*log2(q(x,y)/r(x));
%             end
%         end
%     end
%     C
%     Delta_C=abs(C-C_old);
% end



% iter = 100; % number of iterations
% Nx=size(p,2); % size of input alphabet
% Ny=size(p,1); % size of output alphabet
% 
% %% initialization
% 
% r = rand(1,Nx); 
% r = r/sum(r); % initial input distribution
% 
% q = zeros(Nx,Ny); % q(x,y) 
% 
% %% BAA
% 
% for it=1:iter
% 
% for y=1:Ny % updating q
%     p_sum=0;
%     for x_prime=1:Nx
%         p_sum = p_sum+r(x_prime)*p(y,x_prime);
%     end
%     for x=1:Nx
%         q(x,y) = r(x)*p(y,x)/p_sum;
%     end
% end
% 
% r=ones(1,Nx);
% for x=1:Nx % updating r
%     for y=1:Ny
%         r(x) = r(x)*(q(x,y)^p(y,x));
%     end
% end
% r = r/sum(r);
% 
% end
% 
% C=0;
% for x=1:Nx
%     for y=1:Ny
%         if q(x,y)~=0
%             C = C + r(x)*p(y,x)*log2(q(x,y)/r(x));
%         end
%     end
% end

