%INVERSI BPT
%KELOMPOK 1 TOMOGRAPHY
%DOWNHOLE SEISMIC
%NO BENDING RAYPATH MODEL (2 Layers)
%% Forward Modelling 2 lapis homogen isotrop
clear;close;clc
%Parameter
v1 = 1000,v2 =2000 %kecepatan vp lapisan 
h1 = 50;h2=60; %ketebalan lapisan 1 dan 2
tbl = [h1 h1+h2]; %Tebal Model/Depth batas lapisan
Offset = 0:100; %Lintasan
Oflen = length(Offset);
sh_x  = [ Oflen/4 Oflen/3 Oflen/2 Oflen]; %Shot absis
sh_y  = zeros(size(sh_x));   %Shot ordinat
rec_y = [tbl(2)/4 tbl(2)/3 tbl(2)/2 tbl(2)];  %Receiver depth
rec_x = zeros(size(rec_y));               %Receiver xline
if  h1<=20 || h1>=h2
    error('Titik batas lapisan perlu berada di minimal antara REC 2-3')
end
%LINTASAN RAMBAT SINAR
A = zeros(length(rec_y),length(sh_x)); 
for i=1:length(sh_x)                           
    for j = 1:length(rec_y)
        A(j,i) = sqrt(sh_x(i)^2+rec_y(j)^2);
    end
end
%WAKTU TEMPUH
C = zeros(length(rec_x),length(sh_x));   %Matriks Time Travel
d = rec_y(rec_y>h1);
R = zeros(length(d),length(rec_y));
for i= 1:length(rec_y)
    for j=1:length(d)
        R(j,i) = atand(d(j)/sh_x(i));   %SUDUT RAY TEMBUS BATAS
    end
end
F = zeros(size(R));       %Panjang Ray pada Lapisan 1 Sebelum Tembus Batas
for i=1:length(rec_y)
    for j=1:length(d)
        F(j,i) = h1/sind(R(j,i));
    end
end
Z = A(3:4,:);              %Ekstraksi Lintasan Tembus Lapisan 2
ZZ = Z-F;                  %Lintasan hanya di lapisan 2
U  = zeros(2,4);
ZZ = vertcat(U,ZZ);
FF = vertcat(U,F);
ZZZ = zeros(size(A));      %Matriks Waktu tempuh
for i=1:length(sh_x) 
    for j=1:length(rec_y)
        if j <= length(d)
            ZZZ(j,i) = A(j,i)/v1;
        else 
            ZZZ(j,i) = ZZ(j,i)/v2 + FF(j,i)/v1;   %WAKTU TEMPUH SHOT 1-4 di REC 1-4
        end
    end
end

%PLOTIING MODEL
figure(1)
title('FORWARD MODEL BPT')
ylabel('Depth(m)')
xlabel('Offset(m)')
xlim([0 max(Offset)+10])
ylim([0 max(rec_y)+30])
% plot(h1,0,'-k');hold on
yline(h1,'-k');hold on
xp =[0 0 max(Offset)+10 max(Offset)+10];
yp =[0 max(rec_y)+30 max(rec_y)+30 0];
plot(xp,yp)
Lx1 =[0 0 max(Offset)+10 max(Offset)+10];
Ly1 =[0 h1 h1 0];
Ly2 =[h1 h2+h1+30 h2+h1+30 h1];
plot(Lx1,Ly1);               %LAPISAN 1
fill(Lx1,Ly1,'y');          
hold on
plot(Lx1,Ly2);               %LAPISAN 2
fill(Lx1,Ly2,'g');
hold on
REC = [rec_x' rec_y'];
SH  = [sh_x' sh_y'];
for i=1:length(REC)
    for j = 1:length(SH)
        plot(SH(j,:),REC(i,:),'Marker','V')
        hold on
    end
end
set(gca,'ydir','reverse')
fprintf('AKSES TIME TRAVEL dengan ketik ZZZ, kolom n \nberati shot ke-n, baris n berarti rec ke-n \nLintasan Ray Path ketik A')


% INVERSION BPT
%SLOWNESS RAYPATH (RATA-RATA)

PP = zeros(size(A));

for p=1:length(sh_x)
    for l =1:length(rec_y)
        PP(l,p) = ZZZ(l,p)/A(l,p);     %SLOWNESS TIAP SH(n)-REC(1-4) (RATA2)
    end
end

%CONSTRAIN BATAS
P = zeros(length(sh_x)-1,length(rec_y));
zz1 = sum(ZZZ(1,:))/4;           %WAKTU RATA2
zz2 = sum(ZZZ(2,:))/4;
zz3 = sum(ZZZ(3,:))/4;
zz4 = sum(ZZZ(4,:))/4;
ZX = [zz2-zz1 zz3-zz2 zz4-zz3];

%SLOWNESS LAPISAN KEDUA
SL = zeros(size(PP));
Guest = 5;                           %Tebakan Awal (Increment Constraint Awal dalam permil
S_One = ZZZ(1,1)/A(1,1);              %First Constraints Slowness (SH1-Rec1)
S_Guest= S_One+Guest/1000;

if S_Guest<S_One
    error('S_Guest kurang besar, mungkin v2 akan lebih kecil dari sebenarnya')
elseif S_Guest>S_One+0.005
    error('S_Guest terlalu besar, mungkin v2 akan lebih besar dari sebenarnya')
end
LL = 5;                               %Iterasi (n-1)
dt = 1; 
k=0;
fprintf('\nAKSES SLOWNESS LAPISAN YANG DIINVERSI, KETIK ZL')
if max(ZX) == ZX(1)       %ZONA POTENSI BATAS LAPISAN
    ZY = [rec(1) rec_y(2)];border = 1;ZL = zeros(length(rec_y)-border,2);h0b = max(ZY);
    fprintf('\nZONA POTENSI BATAS DIANTARA REC 1 - REC 2');
    for j=1:2
        for i=1:length(rec_y)-border
            ZL(i,j) = abs(ZZZ(i,j)-S_Guest*A(i+border,j));
        end
    end        
elseif max(ZX) == ZX(2)          %AFTER OR ON rec(border)
    border = 2;
    fprintf('\nZONA POTENSI BATAS DIANTARA REC 2 - REC 3');
    ZY = [rec_y(2) rec_y(3)];ZL = zeros(length(rec_y)-border,2);h0b = max(ZY);
    for j=1:2
        for i=1:length(rec_y)-border
            ZL(i,j) = abs(ZZZ(i,j)-S_Guest*A(i+border,j));
        end
    end        
elseif max(ZX) == ZX(3)
    border = 3;
    fprintf('\nZONA POTENSI BATAS DIANTARA REC 3 - REC 4');
    ZY = [rec_y(3) rec_y(4)];ZL = zeros(length(rec_y)-border,2);h0b = max(ZY);
  for j=1:2
        for i=1:length(rec_y)-border
            ZL(i,j) = abs(ZZZ(i,j)-S_Guest*A(i+border,j));
        end
  end        
end

ZL_Mean = mean(ZL,2);
ZL_Fix  =repmat(ZL,1,2)*S_One*1.4;
LP = PP(1:border,:);
GRID = vertcat(LP,ZL_Fix);
VELO = 1./GRID;
VELO(3,:) = VELO(4,:);
%PLOTTING MODEL TOMO


figure(2)
First_Hor = [0 sh_x max(xp)];
First_Depth = [0 rec_y max(yp)];
a = VELO(1,:);n = VELO(size(VELO,1),:);
Z  = vertcat(a, VELO);Z  = vertcat(Z,n);cc = Z(:,1);
Z  = horzcat(cc,Z,cc);
[X,Y]= meshgrid(First_Hor,First_Depth);
pcolor(X,Y,Z);hold on
colorbar;
set(gca,'Ydir','reverse');
title('Model Kecepatan');
xlabel('Offset Shot(m)');
ylabel('Depth Receiver(m)');
plot(rec_x,rec_y,'Vk','MarkerSize',15)
hold on














