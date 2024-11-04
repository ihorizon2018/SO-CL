%  Structural Truss Problems Benchmark Suite  
%
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
% The truss analyser function is taken from is slightly modified to be used in this code 
% The initial parameters that you need are:
%__________________________________________
% problem_type = 1 to 5, 
%          problem_type -
%            1: 200-truss static optimization
%            2: 10-truss static optimization
%            3: 25-truss static optimization
%            4: 72-truss static optimization
%            5: 120-truss static optimization
% dim = number of your variables
% Max_iteration = maximum number of iterations
% OPTIM_DATA.POP = number of search agents
% OPTIM_DATA.RT = number of running times
% parameter = the control parameters of SO-CL
% OPTIM_DATA.ansname = the filename of optimization results 

% To run SO-CL: main_truss
%______________________________________________________________________________________________


clear all
close all
clc
% Select the problem you want to solve
problem_type = 1;
% Choose the filename for the results
OPTIM_DATA.ansname = 'ansname2';
algorithmType = 'SO_CL';

%% Select your SO-CL parameters
parameter.C1 = 0.5;
parameter.C2 = 0.05;
parameter.C3 = 2;
parameter.Threshold1 = 0.25;
parameter.Thresold2 = 1;
parameter.pv = 0.9;
parameter.Pl = 0.6;

% Set population size, number of runs, and analysis count
OPTIM_DATA.POP = 50;               % Population size
OPTIM_DATA.RT = 2;                 % Runtime
OPTIM_DATA.MXFE_D = 40000/2;         % maxFEs = OPTIM_DATA.MXFE_D * D (FE analysis count / population size)
OPTIM_DATA.name = problem_type;

optimization_solver


save(myfilename,'final');
figure;
plot(mean(curveitSO(:,1:2:end)));hold on
plot(mean(curveit));
legend('SO','SO-CL')
