% This function analysis the design
% The truss analyser function is taken from this link: http://www.mathworks.com/matlabcentral/fileexchange/14313-truss-analysis
% and is slightly modified to be used in this code
function Objection_function=FUN120(Individual)
Truss_Trussdata=TrussTrussdata_120bar;
Truss_Trussdata.A(1:12,:)=Individual(1);
Truss_Trussdata.A(13:24,:)=Individual(2);
Truss_Trussdata.A(25:36,:)=Individual(3);
Truss_Trussdata.A(37:60,:)=Individual(4);
Truss_Trussdata.A(61:84,:)=Individual(5);
Truss_Trussdata.A(85:96,:)=Individual(6);
Truss_Trussdata.A(97:120,:)=Individual(7);

dnsay=size(Truss_Trussdata.dnkoor,1);%单元数目49
elsay=size(Truss_Trussdata.eldn,1);%杆件120


topser=dnsay*3;%自由度数3*49
yer=zeros(topser,1);%生成自由度的零向量
%rijitli刚度矩阵elboy杆件长度
[rijitlik,elboy]=stiffness_3D_truss(topser,elsay,Truss_Trussdata.eldn,Truss_Trussdata.dnkoor,Truss_Trussdata.E,Truss_Trussdata.A);
%tutser约束的自由度号
tutser=find(Truss_Trussdata.MesKos'==1);
%FF外荷载F把FF变为一列
FF=Truss_Trussdata.yuk';
F=FF(:);
%用刚度矩阵算位移yer支座位置为零
ser=setdiff([1:topser]',[tutser]);
    U=inv(rijitlik(ser,ser))*F(ser);
    yer=zeros(topser,1);
    yer(ser)=U;
%gerilmeler应力
[gerilmeler]=stresses_3D_truss(elsay,Truss_Trussdata.eldn,Truss_Trussdata.dnkoor,Truss_Trussdata.E,yer,elboy);


%% Constraints
% Stress Constraint 
maxCEK=0.6*Truss_Trussdata.fy;
alfa=0.4993;
beta=0.6777;
r_i=alfa*(Truss_Trussdata.A).^beta;
k=1;
LAMDA_i=k*elboy'./r_i;
Cc=sqrt(2*pi^2*Truss_Trussdata.E/Truss_Trussdata.fy);

maxBAS=zeros(size(gerilmeler,1),1);

for K=1:size(gerilmeler,1)
    if gerilmeler(K,1)<0
        if LAMDA_i(K) < Cc
            maxBAS(K)=((1-(LAMDA_i(K))^2/(2*Cc^2))*Truss_Trussdata.fy)/(5/3+3*LAMDA_i(K)/(8*Cc)-LAMDA_i(K)^3/(8*Cc^3));
        else
            maxBAS(K)=12*pi^2*Truss_Trussdata.E/(23*(LAMDA_i(K))^2);
        end

        g1(K,1)=(abs(gerilmeler(K,1)))/maxBAS(K)-1;%Basınç
    else
        g1(K,1)=(abs(gerilmeler(K,1)))/maxCEK-1;%Çekme
    end
end
%gg(1)应力超限gg(2)位移超限
gg(1)=sum(g1.*(g1>0));

% Displacement constraint
g2=abs(yer)/Truss_Trussdata.maxYER-1;
gg(2)=10*sum(g2.*(g2>0));

%% Object Function
%ObjVa质量
ObjVal=0;
for i=1:size(Truss_Trussdata.eldn,1)
ObjVal=ObjVal+elboy(i)*Truss_Trussdata.A(i)*Truss_Trussdata.GS;
end

%% Penalized Obj. Func.
%Z罚函数
PEN=10^4;
Z=ObjVal;
for k=1:length(gg)
     Z=Z+ PEN*gg(k);
end
Objection_function=Z;








