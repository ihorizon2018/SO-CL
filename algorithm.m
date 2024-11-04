clear all 
close all
rng("shuffle")
% clc

pop=50; % Number of search agents 种群数量
Max_iteration=500; % Maximum numbef of iterations 设定最大迭代次数

% sphere rastrigin rosenbrock ackley griewank schwefel
%sphere F1 rosenbrock F5 schwefel? F8 rastrigin F9 Ackley F10 griewank F11
% Schwefel_2 f2
% Function_name='rastrigin'; 
% Load details of the selected benchmark function
% [lb,ub,dim,fobj]=My_Functions_details(Function_name);  %设定边界以及优化函数
% figure
% func_plot(Function_name)

Function_name=1; 
fobj = @(x) cec17_func(x',Function_name);
dim=10; %Dimension,[2,10,30,50,100]
lb=-100;%lower boound
ub=100;%upper bound
fobj(ones(1,dim),1)

[Best_posSO,Best_scoreSO,SO_curve] = SO(pop,2*Max_iteration,lb,ub,dim,fobj);

parameter.C1 = 0.5;
parameter.C2 = 0.05;
parameter.C3 = 2;
parameter.Threshold1 = 0.25;
parameter.Thresold2 = 0.6;
parameter.pv = 0.8;
parameter.Pl = 0.5;
[Best_posESO,Best_scoreESO,ESO_curve1]=SO_CL(pop,Max_iteration,lb,ub,dim,fobj,parameter); %开始优化


%% 结果对比


morandi_colors = [
    235,36,38;
    249, 250, 20;
    128, 203, 88;
    39, 150, 235;
    61, 38, 168;
    066, 062, 060;]/255;% 深蓝色

figure('Position',[269   240   500   300])
semilogy(SO_curve(1:2:end),'Color',morandi_colors(5, :),'linewidth',1)
hold on
semilogy(ESO_curve1,'Color',morandi_colors(1, :),'linewidth',1,LineStyle='--')


switch Function_name
    case 'Schwefel_2'
        Function_name='Schwefel function';
    case 'sphere'
        Function_name='Sphere function';
    case 'rastrigin'
        Function_name='Rastrigin function';
    case 'rosenbrock'
        Function_name='Rosenbrock function';
    case 'ackley'
        Function_name='Ackley function';
    case 'griewank'
        Function_name='Griewank function';
    case 'schwefel'
        Function_name='Schwefel function';
    % otherwise
    %     error("wrong input")        
end

xlabel('Iteration');
ylabel('Objective function');
axis normal
% 设置 X 轴为线性等间隔（如果需要，也可以不设置）
xlim([0 Max_iteration]);
ylim([ESO_curve1(end) SO_curve(1)]);
xticks(linspace(0, Max_iteration, 6));  % X 轴等间隔

% 设置 Y 轴的对数刻度等间隔

grid on
box on
legend('SO','MSO')
fontsize=14;
fontname='Arial';
set(gca, 'FontName', fontname, 'FontSize', fontsize); % 修改坐标轴的字体和字号
title(Function_name, 'FontSize', fontsize+2)
set(findall(gca, 'Type', 'text'), 'FontName', fontname, 'FontSize', fontsize); % 修改图中所有文本对象的字体和字号

print(Function_name, '-dpng', '-r400');
print(Function_name, '-djpeg', '-r400');




save([Function_name, '3']);


