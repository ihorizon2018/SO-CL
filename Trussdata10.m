
function D=Trussdata10
Coord=360*[2 1 0;2 0 0;1 1 0;1 0 0;0 1 0;0 0 0]; %节点编号
Con=[5 3;1 3;6 4;4 2;3 4;1 2;5 4;6 3;3 2;4 1];%单元编号
Re=[0 0 1;0 0 1;0 0 1;0 0 1;1 1 1;1 1 1];%节点约束，顺序对应coord,1约束0不约束
Load=zeros(size(Coord));Load(2,:)=[0 -1e5 0];Load(4,:)=[0 -1e5 0];%荷载case 1
E=ones(1,size(Con,1))*1e7;%ea
A=ones(1,10);
% Available sections%可使用界面，对应41种
AV=[1.62, 1.8, 1.99, 2.13, 2.38, 2.62, 2.88, 2.93, 3.09, 3.13, 3.38, 3.47, 3.55, 3.63, 3.84,...
        3.87, 3.88, 4.18, 4.22, 4.49, 4.59, 4.80, 4.97, 5.12, 5.94, 7.22, 7.97, 11.5, 13.50,...
        13.90, 14.2, 15.5, 16.0, 16.9, 18.8, 19.9, 22.0, 22.9, 28.5, 30.0, 33.5];%in^2
%Allowable Stress（容许应力）
TM=25000;%psi
%Allowable Displacement（容许位移）
DM=2;%inch
%WEIGHT PER UNIT of VOLUME（每个单元体积）
RO=.1;%lb/in^3
LB=ones(1,10)*0.1;
UB=ones(1,10)*35;
D=struct('Coord',Coord','Con',Con','Re',Re','Load',Load','E',E','A',A','AV',AV','TM',TM','DM',DM','RO',RO','LB',LB,'UB',UB);%把所有变量都保存到D
%中并返回其变量值