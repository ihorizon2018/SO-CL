%% MAIN CODE
switch problem_type
    case 1
        OPTIM_DATA.ObjectiveFunction = @(Individual) FUN200(Individual);
    case 2
        OPTIM_DATA.ObjectiveFunction = @(Individual) FUN10(Individual);
    case 3
        OPTIM_DATA.ObjectiveFunction = @(Individual) FUN25(Individual);
    case 4
        OPTIM_DATA.ObjectiveFunction = @(Individual) FUN72(Individual);
    case 5
        OPTIM_DATA.ObjectiveFunction = @(Individual) FUN120bar(Individual);
    otherwise
        error("wrong problem_type, please choose 1 to 5");
end
[OPTIM_DATA.LB,OPTIM_DATA.UB]=problem_bounds(problem_type);
%% PARAMETERS
ub=OPTIM_DATA.UB;
lb=OPTIM_DATA.LB;
RT=OPTIM_DATA.RT;
N=OPTIM_DATA.POP; % Population size N
maxT=OPTIM_DATA.MXFE_D/OPTIM_DATA.POP;
fobj=OPTIM_DATA.ObjectiveFunction; % Objective function
dim=size(lb,2); % Variable dimension
myfilename=OPTIM_DATA.ansname;

%% ALGORITHM
% Initialization
% Example assignment (within the loop)
run=1;

min_iter = 0;     % Set the initial value of min_iter as needed
fval = 0;         % Set fval according to the actual optimization problem
cruve_FES = [];   % Initialize as empty, assuming it will be assigned in the loop
curveit = [];     % Initialize curveit for each loop result
bestcruve=[];
algorithmType='SO_CL';
while(run<=RT)


[Xfood(run,:), fval(run), curveit(run,:), cruve_FES(run,:)] = SO_CL(N, maxT, lb(1), ub(1), dim, fobj, parameter);
% [Xfood(run,:), fval(run), curveit(run,:), cruve_FES(run,:)] = SO_CL2(N, maxT, lb(1), ub(1), dim, fobj, parameter);
[XfoodSO(run,:), fvalSO(run), curveitSO(run,:)] = SO(N, maxT*2, lb(1), ub(1), dim, fobj);


min_iter= find(curveit <= round(fval(run), 2), 1, 'first');
if isempty(min_iter)
    min_iter = maxT;  % Assign a default value if empty
end

index(run) =min_iter;
F(run)=fval(run); % Minimum value in this run
MinT(run)=index(run)*N; % Minimum heuristic loop count in this run
X(run,:)=Xfood(run,:); % Size optimization result in this run
if isempty(cruve_FES(run,:))
    bestcruve(run,:)=cruve_FES(run,:);
end
bestcruveit(run,:)=curveit(run,:);
run=run+1;
end
%%%%%%%%% Output useful results
final.Number_of_analyses_per_calculation=MinT; % Minimum analysis count
final.Minimum_number_of_iterations=index; % Minimum iteration count
final.Size_of_each_section=X; % X
final.Quality_results_per_calculation=F; % F
final.OPTIM_DATA=OPTIM_DATA; % Parameters
if isempty(bestcruve)
    final.FEScruve=bestcruve; % FES iteration curve
end

final.cruve=bestcruveit; % FES iteration curve


