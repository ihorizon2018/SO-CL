function [LB,UB]=problem_bounds(problem_type)
% LB : Lower Bounds
% UB : Upper Bounds
% Select Problem
% 1 : 17  Bar Planar Truss Problem
% 2 : 18  Bar Planar Truss Problem
% 3 : 200 Bar Planar Truss Problem
% 4 : 25  Bar Space Truss Problem
% 5 : 72  Bar Space Truss Problem
% 6 : 120 Bar Space Truss Problem
% 7 : 10  Bar Truss Problem with Frequency Constraint
% 8 : 37  Bar Truss Problem with Frequency Constraint
% 9 : 52  Bar Truss Problem with Frequency Constraint

switch problem_type

    case 1
        LB=ones(1,29)*0.1;
        UB=ones(1,29)*20;
    case 2
        LB=ones(1,10)*0.1;
        UB=ones(1,10)*35;
    case 3
        LB=ones(1,8)*0.01;
        UB=ones(1,8)*3.4;
    case 4
        LB=ones(1,16)*0.1;
        UB=ones(1,16)*20;
    case 5
        LB=ones(1,7)*0.775;
        UB=ones(1,7)*20;
    otherwise
        error("wrong problem_type, please choose 1 to 5");


end