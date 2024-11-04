%  Source codes demo version 2.0                                                                      
%                                                                                                     
%  Developed in MATLAB R2024b                                                   %                                                               
%                                                                                                     
%  Author and programmer: Pinghe Ni, Xiaoyu Sua, Jinlong Fu*                                                          
%                                                                                                     
%         e-Mail: sxy123456@emails.bjut.edu.cn                                                             
%                 pinghe.ni@connect.polyu.hk                                                                                                 
%_______________________________________________________________________________________________                                                                                                  
function [Xfood, fval,gbest_t,cruve,count,num] = SO_CL2(N,T,lb,ub,dim,fobj,parameter)
global tweight;
tweight = 1; % 设置 t 的初始值
%% initial
count=1;
count_lens=1;%测试lens利用率
count_criss1=1;
numcriss1=[0];
count_criss=1;
vec_flag=[1,-1];
% Threshold=0.25;
% Thresold2= 0.6;
% C1=0.5;%eq.(14)
% C2=0.05;%eq.(15)

C1=parameter.C1;
C2=parameter.C2;
C3=parameter.C3;
Threshold=parameter.Threshold1;
Thresold2=parameter.Thresold2;
pv=parameter.pv;
Pl=parameter.Pl;


C3=zeros(1,T);
for l= 1:T
    C3(1,l)=2-2*sin(0.5*pi*(l/T)^4);%eq.(16)
end
ub=ub(1);
lb=lb(1);
ub1 = ub.*ones(1,dim);
lb1 = lb.*ones(1,dim);
if(max(size(ub)) == 1)
    ub2 = ub.*ones(1,dim);
    lb2 = lb.*ones(1,dim);
end

%%
X=lb+rand(N,dim)*(ub(1)-lb(1));%eq.(1)
fitness=zeros(1,N);
for i=1:N
    fitness(i)=feval(fobj,X(i,:));
end
Trajectories=zeros(N,T);
position_history=zeros(N,T,dim);
fitness_history=zeros(N,T);


[GYbest, gbest] = min(fitness);
Xfood = X(gbest,:);
%% Diving the swarm into two equal groups males and females
Nm=round(N/2);%eq.(2&3)
Nf=N-Nm;
Xm=X(1:Nm,:);
Xf=X(Nm+1:N,:);
fitness_m=fitness(1:Nm);
fitness_f=fitness(Nm+1:N);
[fitnessBest_m, gbest1] = min(fitness_m);
BestIndex1=gbest1;
Xbest_m = Xm(gbest1,:);
[fitnessBest_f, gbest2] = min(fitness_f);
BestIndex2=gbest2;
Xbest_f = Xf(gbest2,:);
%% % 计算惯性权重因子omega
t = 1:T; % 迭代时间
omega_max = 1.0; % 最大惯性权重
omega_min = 0.2; % 最小惯性权重
% omega = (omega_max - omega_min) .* (1 - t / T) + omega_min;
omega=omega_max.*ones(1,T);

% omega=(exp(2 * (1 - t / T)) - exp(-2 * (1 - t / T))) ./ ...
%     (exp(2 * (1 - t / T)) + exp(-2 * (1 - t / T)));


delta_max = 2; % 根据具体情况设定
% delta_min = 1; % 根据具体情况设定
% k = -1.5*(-1-2*(t/T).^2);% scaling factor eq.(19)
% % k = delta_max/2 .* ((delta_max - delta_min) - 2 .* ((t / T) .^ 0.5));
% delta_max = 1.2;
k = delta_max .* (1 - 2 .* (t / T) .^ 0.01);
% k = 2-2.*(t/T);
% k = -1*ones(1,T);% scaling factor eq.(19)

%% Main loop
for t = 1:T
    %% the principle of convex lens imaging
    tweight=t/T;
    % Pl=rand;
    for i=1:dim
        r3=rand;
        if r3<Pl
            len(i)=(ub + lb)./2 + (ub + lb)./(2*k(t)) - Xbest_m(i)./k(t);
        else
            len(i)=Xbest_m(i);
        end
    end
    % Temp = (ub + lb)./2 + (ub + lb)./(2*k) - Xbest_m./k; %eq.(20)
    Temp=len;
    for c = 1:dim
        if(Temp(1,c)>ub)
            Temp(1,c) =ub1(c);
        end
        if(Temp(1,c)<lb)
            Temp(1,c)  =lb1(c);
        end
    end
    fitTemp = fobj(Temp);
    if(fitTemp<fitnessBest_m)
        fitnessBest_m=fitTemp ;
        Xbest_m = Temp;
        Xm(BestIndex1,:) = Temp;

        count_lens=count_lens+1;%测试lens利用率
        num_lens(count_lens)=t;
    end



    for i=1:dim
        r3=rand;
        if r3<Pl
            len(i)=(ub + lb)./2 + (ub + lb)./(2*k(t)) - Xbest_f(i)./k(t);
        else
            len(i)=Xbest_f(i);
        end
    end
    Temp=len;
    for c = 1:dim
        if(Temp(1,c)>ub)
            Temp(1,c) =ub1(c);
        end
        if(Temp(1,c)<lb)
            Temp(1,c)  =lb1(c);
        end
    end
    fitTemp = fobj(Temp);
    if(fitTemp<fitnessBest_f)
        fitnessBest_f=fitTemp ;
        Xbest_f = Temp;
        Xf(BestIndex2,:) = Temp;

        count_lens=count_lens+1;%测试lens利用率
        num_lens(count_lens)=t;
    end

    Temp=exp(-((t)/T));  %eq.(4)
    Q=C1*exp(((t-T)/(T)));%eq.(5)

    if Q>1        Q=1;    end
    %% Exploration Phase (No Food)
    if Q<Threshold
        for i=1:Nm
            for j=1:1:dim
                rand_leader_index = floor(Nm*rand()+1);
                X_randm = Xm(rand_leader_index, :);
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                Am=exp(-fitness_m(rand_leader_index)/(fitness_m(i)+eps));
                Xnewm(i,j)=X_randm(j)+Flag*C2*Am*((ub-lb)*rand+lb);%eq.(5)
            end
        end
        for i=1:Nf
            for j=1:1:dim
                rand_leader_index = floor(Nf*rand()+1);
                X_randf = Xf(rand_leader_index, :);
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                Af=exp(-fitness_f(rand_leader_index)/(fitness_f(i)+eps));
                Xnewf(i,j)=X_randf(j)+Flag*C2*Af*((ub-lb)*rand+lb);%eq.(6)
            end
        end
        %% Exploitation Phase (Food Exists)
    else
        if Temp>Thresold2  %hot
            for i=1:Nm
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                for j=1:1:dim
                    Xnewm(i,j)=omega(t).*Xfood(j)+C3(1,t)*Flag*Temp*rand*(Xfood(j)-Xm(i,j));%eq.(7)
                end
            end
            for i=1:Nf
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                for j=1:1:dim
                    Xnewf(i,j)=omega(t).*Xfood(j)+Flag*C3(1,t)*Temp*rand*(Xfood(j)-Xf(i,j));%eq.(7)
                end
            end
        else %cold
            if rand>0.6 %fight
                for i=1:Nm
                    for j=1:1:dim
                        FM=exp(-(fitnessBest_f)/(fitness_m(i)+eps));
                        Xnewm(i,j)=Xm(i,j) +C3(1,t)*FM*rand*(omega(t).*Q*Xbest_f(j)-Xm(i,j));%eq.(8)

                    end
                end
                for i=1:Nf
                    for j=1:1:dim
                        FF=exp(-(fitnessBest_m)/(fitness_f(i)+eps));
                        Xnewf(i,j)=Xf(i,j)+C3(1,t)*FF*rand*(omega(t).*Q*Xbest_m(j)-Xf(i,j));%eq.(9)
                    end
                end
            else%mating
                for i=1:Nm
                    for j=1:1:dim
                        Mm=exp(-fitness_f(i)/(fitness_m(i)+eps));
                        Xnewm(i,j)=Xm(i,j) +C3(1,t)*rand*Mm*(Q*Xf(i,j)-Xm(i,j));%eq.(10)
                    end
                end
                for i=1:Nf
                    for j=1:1:dim
                        Mf=exp(-fitness_m(i)/(fitness_f(i)+eps));
                        Xnewf(i,j)=Xf(i,j) +C3(1,t)*rand*Mf*(Q*Xm(i,j)-Xf(i,j));%eq.(11)
                    end
                end
            end
        end
    end
    %% Return back the search agents that go beyond the boundaries of the search space
    for j=1:Nm
        Flag4ub=Xnewm(j,:)>ub;
        Flag4lb=Xnewm(j,:)<lb;
        Xnewm(j,:)=(Xnewm(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,Xnewm(j,:));
        if y<fitness_m(j)
            fitness_m(j)=y;
            Xm(j,:)= Xnewm(j,:);
        end
        cruve(count)=GYbest;
        count=count+1;
    end

    %% Return back the search agents that go beyond the boundaries of the search space
    for j=1:Nf
        Flag4ub=Xnewf(j,:)>ub;
        Flag4lb=Xnewf(j,:)<lb;
        Xnewf(j,:)=(Xnewf(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,Xnewf(j,:));
        if y<fitness_f(j)
            fitness_f(j)=y;
            Xf(j,:)= Xnewf(j,:);
        end
        cruve(count)=GYbest;
        count=count+1;
    end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%十字交叉搜寻%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%%%%%%%%%%%%%%%% 平行交叉
    % 雄性
    for j = 1:10
        % 改进2.1：★★横向交叉策略★★
        if (mod(j,2) == 1)
            for i = 1 : dim
                r = rand();
                c = 2*rand()-1;
                xm_criss( j, i ) = r * Xm( j, i )+(1-r)*Xm( j+1, i )+c*(Xm( j, i )-Xm( j+1, i ));
            end
        else
            for i = 1 : dim
                r = rand();
                c = 2*rand()-1;
                xm_criss(j, i ) = r * Xm( j-1, i )+(1-r)*Xm( j, i )+c*(Xm( j-1, i )-Xm( j, i ));
            end
        end
        Flag4ub=xm_criss(j,:)>ub1;
        Flag4lb=xm_criss(j,:)<lb1;
        xm_criss(j,:)=(xm_criss(j,:).*(~(Flag4ub+Flag4lb)))+ub1.*Flag4ub+lb1.*Flag4lb;
        y = feval(fobj,xm_criss(j,:));
        if y<fitness_m(j)
            fitness_m(j)=y;
            Xm(j,:)= xm_criss(j,:);
            count_criss=count_criss+1;%测试criss利用率
            numcriss(count_criss)=count;
        end
        cruve(count)=GYbest;
        count=count+1;
    end 

    % 雌性
    for j = 1:10
        % 改进2.1：★★横向交叉策略★★
        if (mod(j,2) == 1)
            for i = 1 : dim
                r = rand();
                c = 2*rand()-1;
                xf_criss( j, i ) = r * Xf( j, i )+(1-r)*Xf( j+1, i )+c*(Xf( j, i )-Xf( j+1, i ));
            end
        else
            for i = 1 : dim
                r = rand();
                c = 2*rand()-1;
                xf_criss(j, i ) = r * Xf( j-1, i )+(1-r)*Xf( j, i )+c*(Xf( j-1, i )-Xf( j, i ));
            end
        end
        Flag4ub=xf_criss(j,:)>ub1;
        Flag4lb=xf_criss(j,:)<lb1;
        xf_criss(j,:)=(xf_criss(j,:).*(~(Flag4ub+Flag4lb)))+ub1.*Flag4ub+lb1.*Flag4lb;
        y = feval(fobj,xf_criss(j,:));
        if y<fitness_f(j)
            fitness_f(j)=y;
            Xf(j,:)= xf_criss(j,:);
            count_criss=count_criss+1;%测试criss利用率
            numcriss(count_criss)=count;
        end
        cruve(count)=GYbest;
        count=count+1;
    end
    %%%%%%%%%%%%%%%%%%垂直交叉
    for j = 1:10
        % 改进2.2：★★纵向交叉策略★★
        for i = 1 : dim
            d1 = randperm(dim,1);
            d2 = randperm(dim,1);
            r1 = rand();
            x_cross(j, i) = r1*Xm( j, d1 ) + (1-r1)*Xm(j, d2 );
        end
        Flag4ub=x_cross(j,:)>ub1;
        Flag4lb=x_cross(j,:)<lb1;
        x_cross(j,:)=(x_cross(j,:).*(~(Flag4ub+Flag4lb)))+ub1.*Flag4ub+lb1.*Flag4lb;
        y = feval(fobj,x_cross(j,:));
        if y<fitness_m(j)
            fitness_m(j)=y;
            Xm(j,:)= x_cross(j,:);
            count_criss=count_criss+1;%测试criss利用率
            numcriss(count_criss)=count;
        end
        cruve(count)=GYbest;
        count=count+1;
    end
    for j = 1:10
        % 改进2.2：★★纵向交叉策略★★
        for i = 1 : dim
            d1 = randperm(dim,1);
            d2 = randperm(dim,1);
            r1 = rand();
            x_cross(j, i) = r1*Xf( j, d1 ) + (1-r1)*Xf(j, d2 );            
        end
        Flag4ub=x_cross(j,:)>ub1;
            Flag4lb=x_cross(j,:)<lb1;
            x_cross(j,:)=(x_cross(j,:).*(~(Flag4ub+Flag4lb)))+ub1.*Flag4ub+lb1.*Flag4lb;
            y = feval(fobj,x_cross(j,:));
            if y<fitness_f(j)
                fitness_f(j)=y;
                Xf(j,:)= x_cross(j,:);
                count_criss=count_criss+1;%测试criss利用率
                numcriss(count_criss)=count;
            end
            cruve(count)=GYbest;
            count=count+1;
    end

    [Ybest1,gbest1] = min(fitness_m);
    BestIndex1=gbest1;
    [Ybest2,gbest2] = min(fitness_f);
    BestIndex2=gbest2;
    if Ybest1<fitnessBest_m
        Xbest_m = Xm(gbest1,:);
        fitnessBest_m=Ybest1;
    end
    if Ybest2<fitnessBest_f
        Xbest_f = Xf(gbest2,:);
        fitnessBest_f=Ybest2;

    end
    if Ybest1<Ybest2
        gbest_t(t)=min(Ybest1);
    else
        gbest_t(t)=min(Ybest2);

    end
    if fitnessBest_m<fitnessBest_f
        GYbest=fitnessBest_m;
        Xfood=Xbest_m;
    else
        GYbest=fitnessBest_f;
        Xfood=Xbest_f;
    end
    % [t GYbest Xfood]
    [t GYbest]
end
fval = GYbest;
count(1)=count_lens;
count(2)=count_criss;
count(3)=count_criss1;
% num=[num_lens,numcriss];
num=numcriss;

["SO-CL2"]
[count GYbest]
% pause

end





