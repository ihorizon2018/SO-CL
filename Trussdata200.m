function D=Trussdata200
%英寸版
yh=fliplr([0 linspace(360,360+1440,11)]);
Coord=[];
for i=1:11
    if floor(i/2)~=i/2
        Coord=[Coord;[linspace(0,960,5)' yh(i)*ones(5,1)]];
    else
        Coord=[Coord;[linspace(0,960,9)' yh(i)*ones(9,1)]];
    end
end
Coord=[Coord;[240 0];[720 0]];
Coord(1:length(Coord),3)=0;

%  Connectivity单元信息
%%
ele=[];
for i=1:5
    n1=14*i-13;n2=14*i-12;n3=14*i-11;n4=14*i-10;
    n5=14*i-9;n6=14*i-8;n7=14*i-7;n8=14*i-6;
    n9=14*i-5;n10=14*i-4;n11=14*i-3;n12=14*i-2;
    n13=14*i-1;n14=14*i-0;
    n15=n14+1;n16=n15+1;n17=n16+1;n18=n17+1;n19=n18+1;
    elei1=[n1 n2
        n2 n3
        n3 n4
        n4 n5
        n1 n6
        n1 n7
        n2 n7
        n2 n8
        n2 n9
        n3 n9
        n3 n10
        n3 n11
        n4 n11
        n4 n12
        n4 n13
        n5 n13
        n5 n14];
    elei2=[n6 n7
        n7 n8
        n8 n9
        n9 n10
        n10 n11
        n11 n12
        n12 n13
        n13 n14
        n6 n15
        n7 n15
        n7 n16
        n8 n16
        n9 n16
        n9 n17
        n10 n17
        n11 n17
        n11 n18
        n12 n18
        n13 n18
        n13 n19
        n14 n19];
    ele=[ele;elei1;elei2];
end
Con=[ele
    71 72
    72 73
    73 74
    74 75
    71 76
    72 76
    73 76
    73 77
    74 77
    75 77];


%%

% Definition of Degree of freedom (free=0 &  fixed=1); for 2-D trusses the last column is equal to 1
Re=zeros(size(Coord));Re(76:77,:)=[1 1 1;1 1 1 ;];%1约束0不约束
Re(1:length(Coord),3)=1;
% Definition of Nodal loads
loadx=1e3;%N单位牛，对应1000磅，对应千克米每二次方秒，所以用这些单位算出来千克
loady=1e4;%N单位牛，对应1万磅
% loadx=1000;%ip:one  pounds force
% loady=10000;%ip
LoadCase1=zeros(size(Coord));
LoadCase1([1 6  15  20  29 34  43  48  57  62 71],:)=[
    loadx,0,0;
    loadx,0,0;
    loadx,0,0;
    loadx,0,0;
    loadx,0,0;
    loadx,0,0;
    loadx,0,0;
    loadx,0,0;
    loadx,0,0;
    loadx,0,0;
    loadx,0,0;];
for ii=1:55
    loady2(ii,:)=[0,-(loady),0;];
end
LoadCase2=zeros(size(Coord));
LoadCase2([1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 15, ...
    16, 17, 18, 19, 20, 22, 24, 26, 28, 29, 30, ...
    31, 32, 33, 34, 36, 38, 40, 42, 43, 44, 45, 46, ...
    47, 48, 50, 52, 54, 56, 57, 58, 59, 60, 61, 62, 64, ...
    66, 68, 70, 71, 72, 73, 74,75],:)=loady2;
LoadCase3=zeros(size(Coord));
LoadCase3(:,1)=LoadCase1(:,1);
LoadCase3(:,2)=LoadCase2(:,2);
% Definition of Modulus of Elasticity
E=ones(1,size(Con,1))*3e7;%psi
% E=ones(1,size(Con,1))*206e6;%Gpa

% Initialize Section Areas
A=ones(1,length(Con));

%Allowable stress
% TM=2.07e8;%=30 000Ksi;pa=N/m2
TM=1e4;%Ksi
% RO=0.283599305;%lb/in^3
RO=0.283;%lb/in^3

% Groups
Group={{[1:4];[5, 8, 11, 14, 17];[19, 20, 21, 22, 23, 24];[18, 25, 56, 63, 94, 101, 132, 139, 170, 177];
    [26, 29, 32, 35, 38];[6, 7, 9, 10, 12, 13, 15, 16, 27, 28, 30, 31, 33, 34, 36, 37];[39, 40, 41, 42];
    [43, 46, 49, 52, 55];[57, 58, 59, 60, 61, 62];[64, 67, 70, 73, 76],...
    ;[44, 45, 47, 48, 50, 51, 53, 54, 65, 66, 68, 69, 71, 72, 74,75];[77, 78, 79, 80];
    [81, 84, 87, 90, 93];[95, 96, 97, 98, 99, 100];[102, 105, 108, 111, 114];
    [82, 83, 85, 86, 88, 89, 91, 92, 103, 104, 106, 107, 109, 110, 112, 113];[115, 116, 117, 118];
    [119, 122, 125, 128, 131];[133, 134, 135, 136, 137, 138];[140, 143, 146, 149, 152];
    [120, 121, 123, 124, 126, 127, 129, 130, 141, 142, 144, 145, 147, 148, 150, 151];
    [153, 154, 155, 156];[157, 160, 163, 166, 169];[171, 172, 173,174, 175, 176];[178, 181, 184, 187, 190];
    [158, 159, 161, 162, 164, 165, 167, 168, 179, 180, 182, 183, 185, 186, 188,189];[191, 192, 193, 194];
    [195, 197, 198, 200];[196, 199];}};


% Convert to structure array
D=struct('Coord',Coord','Con',Con','Re',Re','LoadCase1',LoadCase1','LoadCase2',LoadCase2','LoadCase3',LoadCase3', ...
    'E',E','A',A','TM',TM','RO',RO','Group',Group);

%% 把塔架画出来
%    h=figure;
% set(gcf,'Position',[0,0,500,1000], 'color','w')
%    box off;axis off;
% for iElem = 1:length(Con(:,1))
%     indexPoint = Con(iElem,:);
% 
%     plot3(Coord(indexPoint,1),Coord(indexPoint,2),Coord(indexPoint,3),'-o','Color','k','MarkerSize',10,...
%     'MarkerFaceColor','#D9FFFF')
%     hold on
% end
% axis off

















