% This function analysis the design
% The truss analyser function is taken from this link: http://www.mathworks.com/matlabcentral/fileexchange/14313-truss-analysis
% and is slightly modified to be used in this code
function GOAL=FUN10(In)%这个是罚函数计算过程
global tweight; % 声明全局变量 t
D=Trussdata10;
D.A=In(1:10)';
% xi2=In(11);
w=size(D.Re);
S=zeros(3*w(2));%刚度矩阵
U=1-D.Re;%节点位移
f=find(U);%寻找自由度f坐标
WE=0;
for i=1:size(D.Con,2)
    H=D.Con(:,i);%单元节点编号
    C=D.Coord(:,H(2))-D.Coord(:,H(1));%单元长度向量
    Le=norm(C);%杆件长度
    T=C/Le;%归一化向量or方向向量cosx sinx
    s=T*T';%两个单元转换矩阵乘积（但这个矩阵s也是四分之一个）
    G=D.E(i)*D.A(i)/Le;%EA/L
    Tj(:,i)=G*T;%力的转置矩阵
    e=[3*H(1)-2:3*H(1),3*H(2)-2:3*H(2)];%整体坐标
    S(e,e)=S(e,e)+G*[s -s;-s s];%刚度的转置矩阵
    WE=WE+Le*D.A(i)*D.RO;%质量
end
U(f)=S(f,f)\D.Load(f );%节点位移，矩阵位移法求解完毕
F=sum(Tj.*(U(:,D.Con(2,:))-U(:,D.Con(1,:))));%整体坐标下的力
% s1=Tj(:,2).*sss(:,2)
% sum(s1)
R=reshape(S*U(:),w);%S*U刚度乘以位移等于杆端弯矩
R(f)=0;%没约束地方没有力（f），实际上只有支座地方有力
TS=(((abs(F'))./D.A)/D.TM)-1;%Tension   %%超出应力的罚函数
US=abs(U')/D.DM-1;%Displacement    %%超出位移的罚函数
PS=sum(TS.*(TS>0));     %%(TS>0)大于零部分为1，小于为0
PD=sum(sum(US.*(US>0)));%同理
GOAL=WE*(1+0.2*PS+1*PD).^(1+1*tweight);
