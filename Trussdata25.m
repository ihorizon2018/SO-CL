function D=Trussdata25
% （data25表示连续的结构截面）
%  Nodal Coordinates（节点：三维空间点）
Coord=[-37.5 0 200;37.5 0 200;-37.5 37.5 100;37.5 37.5 100;37.5 -37.5 100;-37.5 -37.5 100;-100 100 0;100 100 0;100 -100 0;-100 -100 0];

%  Connectivity（连接顺序：从上至下）
Con=[1 2;1 4;2 3;1 5;2 6;2 4;2 5;1 3;1 6;3 6;4 5;3 4;5 6;3 10;6 7;4 9;5 8;4 7;3 8;5 10;6 9;6 10;3 7;4 8;5 9];

% Definition of Degree of freedom (free=0 &  fixed=1); for 2-D trusses the last column is equal to 1
% （定义自由度）（7-10节点为固定端）
Re=zeros(size(Coord));Re(7:10,:)=[1 1 1;1 1 1;1 1 1;1 1 1];

% Definition of Nodal loads, different load cases（定义节点荷载，共有两种工况类型）
LoadCase1=zeros(size(Coord));LoadCase1([1:3,6],:)=1e3*[1,10,-5;0,10,-5;0.5,0,0;0.5,0,0];
LoadCase2=zeros(size(Coord));LoadCase2([1:2],:)=1e3*[0,20,-5;0,-20,-5];

% Define modulus of elasticity（弹性模量）
E=ones(1,size(Con,1))*1e7;

% Initialize the section area（初始界面域）
A=zeros(1,size(Con,1)); 

% Available sections（可挑选截面:连续值）
AV=[.1*[1:26],2.8,3,3.2,3.4,3.2,3];%in^2
% Allowable compression, tension
TMC=[35.092;11.590;17.305;35.092;35.092;6.759;6.959;11.082]*1000;
TMT=40*1000;
% Allowable displacement
DM=0.35;%inch
%WEIGHT PER UNIT VOLUME
RO=.1;%lb/in^3
% Grouping（相同长度单元杆件分为一组，共八组）
Group={{[1];[2;3;4;5];[6;7;8;9];[10;11];[12;13];[14;15;16;17];[18;19;20;21];[22;23;24;25]}};

LB=ones(1,8)*0.01;
UB=ones(1,8)*3.4;
D=struct('Coord',Coord','Con',Con','Re',Re','E',E','A',A','AV',AV','DM',DM','RO',RO','Group',Group,'LB',LB,'UB',UB,'TMC',TMC,'TMT',TMT,'LoadCase1',LoadCase1','LoadCase2',LoadCase2');
