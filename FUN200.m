% This function analysis the design
% The truss analyser function is taken from this link: http://www.mathworks.com/matlabcentral/fileexchange/14313-truss-analysis
% and is slightly modified to be used in this code
function Goal=FUN200(In)
global tweight; % 声明全局变量 t
D=Trussdata200;
for J=1:size(D.Group,1)
    In1(D.Group{J},:)=In(J);
end
D.A=In1;

w=size(D.Re);S=zeros(3*w(2));U=1-D.Re;f=find(U);
WE=0;
for i=1:size(D.Con,2)
   H=D.Con(:,i);%单元节点编号
   C=D.Coord(:,H(2))-D.Coord(:,H(1));%单元长度向量
   Le=norm(C);%杆件长度
   T=C/Le;%归一化向量or方向向量cosx sinx
   s=T*T';%四分之一个转置矩阵
   G=D.E(i)*abs(D.A(i))/Le;%EA/L
   Tj(:,i)=G*T;%力的转置矩阵
   e=[3*H(1)-2:3*H(1),3*H(2)-2:3*H(2)];%整体坐标
   S(e,e)=S(e,e)+G*[s -s;-s s];%刚度的转置矩阵
   WE=WE+Le*abs(D.A(i))*D.RO;
end
%analyse for LoadCase1
U(f)=S(f,f)\D.LoadCase1(f);
F=sum(Tj.*(U(:,D.Con(2,:))-U(:,D.Con(1,:))));%桁架每根杆件内力
R=reshape(S*U(:),w);
R(f)=0;
TS1=(((abs(F'))./D.A)/D.TM)-1;%tanesh

%analyse for LoadCase2
U(f)=S(f,f)\D.LoadCase2(f);
F=sum(Tj.*(U(:,D.Con(2,:))-U(:,D.Con(1,:))));%结果：杆端力
R=reshape(S*U(:),w);R(f)=0;
TS2=(((abs(F'))./D.A)/D.TM)-1;%Tension

%analyse for LoadCase3
U(f)=S(f,f)\D.LoadCase3(f);
F=sum(Tj.*(U(:,D.Con(2,:))-U(:,D.Con(1,:))));
R=reshape(S*U(:),w);R(f)=0;
TS3=(((abs(F'))./D.A)/D.TM)-1;%Tension



TS=max([TS1,TS2,TS3],[],2);
PS=sum(TS.*(TS>0));
GOAL=WE+0.1*WE*PS*(1+1*tweight).^(1+1*tweight);
% Goal=WE+2e3*PS*(1+10*tweight);

