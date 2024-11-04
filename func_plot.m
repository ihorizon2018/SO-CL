% This function draw the benchmark functions
% rosenbrock F5 schwefel F8 rastrigin F9 Ackley F10 griewank F11
function func_plot(Function_name)

[lb,ub,dim,fun]=My_Functions_details(Function_name);

switch Function_name 
    case 'sphere' 
        x=-100:2:100; y=x; %[-100,100]
        
    case 'Schwefel_2' 
        x=-100:2:100; y=x; %[-10,10]
        
    case 'F3' 
        x=-100:2:100; y=x; %[-100,100]
        
    case 'F4' 
        x=-100:2:100; y=x; %[-100,100]
    case 'rosenbrock' 
        x=-200:2:200; y=x; %[-5,5]
    case 'F6' 
        x=-100:2:100; y=x; %[-100,100]
    case 'F7' 
        x=-1:0.03:1;  y=x  %[-1,1]
    case 'schwefel' 
        x=-500:10:500;y=x; %[-500,500]
    case 'rastrigin' 
        x=-5:0.1:5;   y=x; %[-5,5]    
    case 'Ackley' 
        x=-20:0.5:20; y=x;%[-500,500]
    case 'griewank' 
        x=-500:10:500; y=x;%[-0.5,0.5]
    case 'F12' 
        x=-10:0.1:10; y=x;%[-pi,pi]
    case 'F13' 
        x=-5:0.08:5; y=x;%[-3,1]
    case 'F14' 
        x=-100:2:100; y=x;%[-100,100]
    case 'F15' 
        x=-5:0.1:5; y=x;%[-5,5]
    case 'F16' 
        x=-1:0.01:1; y=x;%[-5,5]
    case 'F17' 
        x=-5:0.1:5; y=x;%[-5,5]
    case 'F18' 
        x=-5:0.06:5; y=x;%[-5,5]
    case 'F19' 
        x=-5:0.1:5; y=x;%[-5,5]
    case 'F20' 
        x=-5:0.1:5; y=x;%[-5,5]        
    case 'F21' 
        x=-5:0.1:5; y=x;%[-5,5]
    case 'F22' 
        x=-5:0.1:5; y=x;%[-5,5]     
    case 'F23' 
        x=-5:0.1:5; y=x;%[-5,5]  
end    

    

L=length(x);
f=[];

for i=1:L
    for j=1:L
        if strcmp(Function_name,'F15')==0 && strcmp(Function_name,'F19')==0 && strcmp(Function_name,'F20')==0 && strcmp(Function_name,'F21')==0 && strcmp(Function_name,'F22')==0 && strcmp(Function_name,'F23')==0
            f(i,j)=fun([x(i),y(j)]);
        end
        if strcmp(Function_name,'F15')==1
            f(i,j)=fun([x(i),y(j),0,0]);
        end
        if strcmp(Function_name,'F19')==1
            f(i,j)=fun([x(i),y(j),0]);
        end
        if strcmp(Function_name,'F20')==1
            f(i,j)=fun([x(i),y(j),0,0,0,0]);
        end       
        if strcmp(Function_name,'F21')==1 || strcmp(Function_name,'F22')==1 ||strcmp(Function_name,'F23')==1
            f(i,j)=fun([x(i),y(j),0,0]);
        end          
    end
end

surfc(x,y,f,'LineStyle','none');

end

