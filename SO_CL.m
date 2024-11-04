%  Source codes demo version 1.0                                                                      
%                                                                                                     
%  Developed in MATLAB R2024b                                                   %                                                               
%                                                                                                     
%  Author and programmer: Pinghe Ni, Xiaoyu Sua, Jinlong Fu*                                                          
%                                                                                                     
%         e-Mail: sxy123456@emails.bjut.edu.cn                                                             
%                 pinghe.ni@connect.polyu.hk
%
%  Main paper: A hybrid snake optimizer with crisscross learning strategy       %
%  for constrainedstructural optimization                                       %
                                                                                                    
%_______________________________________________________________________________________________
function [Xfood, fval,curve_it,cruve,count] = SO_CL(N,T,lb,ub,dim,fobj,parameter)
global tweight;
tweight = 1; % 设置 t 的初始值
%% initial
numcriss=1;
count=1;
count_lens=1;%测试lens利用率
  count_SO=1;
count_criss=1;
vec_flag=[1,-1];

C1=parameter.C1;
C2=parameter.C2;
C3=parameter.C3;
Threshold=parameter.Threshold1;
Threshold2=parameter.Thresold2;
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
X=lb+rand(N,dim)*(ub-lb);%eq.(1)
fitness=zeros(1,N);
for i=1:N
    fitness(i)=feval(fobj,X(i,:));
end
Trajectories=zeros(N,T);
position_history=zeros(N,T,dim);
fitness_history=zeros(N,T);
% t1=zeros(1,T);
% t2=zeros(1,T);
% a=0.05;


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

omega=omega_max.*ones(1,T);

delta_max = 2; % 根据具体情况设定

k = delta_max .* (1 - 2 .* (t / T) .^ 0.01);

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
        if Temp>Threshold2  %hot
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
                [~, index]=sort(fitness_m);
                [~, index1]= sort(fitness_f);%排序
                % Xnewm(index(end-3),:)=Xm(index(end-3),:).*(1 + tan(pi*(rand-0.5)));
                Xnewm(index(end-2),:)=Xm(index(end-2),:).*(1 + tan(pi*(rand-0.5)));
                Xnewm(index(end-1),:)=Xm(index(end-1),:).*(1 + tan(pi*(rand-0.5)));
                Xnewm(index(end),:)=Xm(index(end),:).*(1 + tan(pi*(rand-0.5)));
                % Xnewf(index1(end-3),:)=Xf(index1(end-3),:).*(1 + tan(pi*(rand-0.5)));
                Xnewf(index1(end-1),:)=Xf(index1(end-2),:).*(1 + tan(pi*(rand-0.5)));
                Xnewf(index1(end-2),:)=Xf(index1(end-1),:).*(1 + tan(pi*(rand-0.5)));
                Xnewf(index1(end),:)=Xf(index1(end),:).*(1 + tan(pi*(rand-0.5)));
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
            count_SO=count_SO+1;%测试SO利用率
            numSO(count_SO)=count;
        end
        cruve(count)=GYbest;
        count=count+1;
    end
    for j=1:Nf
        Flag4ub=Xnewf(j,:)>ub;
        Flag4lb=Xnewf(j,:)<lb;
        Xnewf(j,:)=(Xnewf(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,Xnewf(j,:));
        if y<fitness_f(j)
            fitness_f(j)=y;
            Xf(j,:)= Xnewf(j,:);
            count_SO=count_SO+1;%测试SO利用率
            numSO(count_SO)=count;
        end
        cruve(count)=GYbest;
        count=count+1;
    end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%交叉基于的综合学习策略%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%crossover-based comprehensive learning (CCL)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLS
    %male
    B=randperm(Nm);
    for j = 1:Nm
        for i=1:dim
            P=rand;
            if(P < pv)
                %% 平行交叉
                r = rand;
                c = 2*rand-1;
                nol=B(j);
                Temp(i) =r*Xm(j,i)+(1-r)*Xm(nol,i)+c*(Xm(j,i)-Xm(nol,i));
            else
                %% 垂直交叉
                ii=randperm(dim,1);
                r4=rand;
                Temp(i)=r4*Xm(j,i)+(1-r4)*Xm(j,ii);
            end
        end
        Temp(Temp>ub) = ub2(Temp>ub);
        Temp(Temp<lb) = lb2(Temp<lb);
        ftemp = fobj(Temp);
        if(ftemp<fitness_m(j))
            fitness_m(j)= ftemp;
            Xm(j,:) = Temp;
            count_criss=count_criss+1;%测试criss利用率
            numcriss(count_criss)=count;
        end
        cruve(count)=GYbest;
        count=count+1;
    end

    %female
    B=randperm(Nf);
    for j = 1:Nf
        for i=1:dim
            P=rand;
            if(P < pv)
                %% 平行交叉
                r = rand;
                c = 2*rand-1;
                nol=B(j);
                Temp(i) =r*Xf(j,i)+(1-r)*Xf(nol,i)+c*(Xf(j,i)-Xf(nol,i));
            else
                %% 垂直交叉
                ii=randperm(dim,1);
                r4=rand;
                Temp(i)=r4*Xf(j,i)+(1-r4)*Xf(j,ii);
            end
        end
        Temp(Temp>ub) = ub2(Temp>ub);
        Temp(Temp<lb) = lb2(Temp<lb);
        ftemp = fobj(Temp);
        if(ftemp<fitness_f(j))
            fitness_f(j)= ftemp;
            Xf(j,:) = Temp;
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
        curve_it(t)=min(Ybest1);
    else
        curve_it(t)=min(Ybest2);

    end
    if fitnessBest_m<fitnessBest_f
        GYbest=fitnessBest_m;
        Xfood=Xbest_m;
    else
        GYbest=fitnessBest_f;
        Xfood=Xbest_f;
    end
    [T GYbest]

end
fval = GYbest;
count(1)=count_lens;
count(2)=count_criss;
 count(3)= count_SO;
% num=[num_lens,numcriss];
num=[numcriss];

["SO-CL 实验组"]
[count GYbest pv]
% pause

end





