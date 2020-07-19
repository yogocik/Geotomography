%GEOTOMOGRAFI ART INVERSION
%KELOMPOK 1 
%FIXED MODEL (4 layers and 6 receivers)
%SHORT INFO---->VARIASI HANYA PADA JUMLAH SHOT,OFFSET,KECEPATAN LAPISAN,DAN
%KETEBALAN (Dengan syarat tebal lapisan pertama tidak boleh kurang dari
%posisi receiver 1 dan jumlah semua tebal lapisan sama dengan Depth)
clear;clc;close all
Depth = 30;                  
Offset = 50;            
R = linspace(5,Depth,6);            % Receiver Offset
n_shot = 40;                        % JUMLAH SHOT
shot = linspace(0,Offset,n_shot);   % Shot Offset
tbl = [5,5,5,15];                  % Tebal per lapisan
w  = 0;                             % Random Noise Ratio
v = [3500,1000,5000,2000];
bounds = [tbl(1) tbl(1)+tbl(2) tbl(1)+tbl(2)+tbl(3) sum(tbl)];     %Batas Lapisan
if length(tbl)~= 4
    error('Jumlah lapisan harus empat')
end
if Depth ~=30
    error('Depth tidak boleh diubah (Depth maximum = 30)')
end
if max(bounds) ~= Depth
    error('Jumlah tebal lapisan tidak sama dengan depth model fix')
end
if min(bounds)<min(R)
    error('Receiver 1 harus berada sebelum/tepat di batas lapisan 1(bounds 1)')
end
if length(v)~=length(tbl)
    error('Jumlah kecepatan tidak sama dengan jumlah lapisan')
end
%RAYPATH
RR = zeros(length(R),length(shot));
for i = 1:length(shot)
    for j = 1:length(R)
        RR(j,i) = ((shot(i))^2+(R(j))^2)^0.5;
    end
end
%TRAVEL TIMES
TT = zeros(size(RR));
fprintf('TRAVEL TIMES MODEL FORWARD (PURE)-->TTA\nRAYPATH LENGTH MODEL FORWARD-->RR\n')
fprintf('TRAVEL TIMES MODEL FORWARD (NOISE)-->TT\n')
for i = 1:length(shot)
    for j = 1:length(R)
        if R(j)<=bounds(1)
            TT(j,i) = RR(j,i)/v(1);
        elseif R(j)<=bounds(2)
            TT(j,i) = ((R(j)-bounds(1))/R(j))*RR(j,i)/v(2) + (bounds(1)/R(j))*RR(j,i)/v(1);
        elseif R(j)<=bounds(3)
            TT(j,i) = (tbl(1)/R(j))*RR(j,i)/v(1) + (tbl(2)/R(j))*RR(j,i)/v(2) + ((R(j)-tbl(2)-tbl(1))/R(j))*RR(j,i)/v(3);
        else
            TT(j,i) = ((tbl(1)/R(j))*RR(j,i))/v(1) + ((tbl(2)/R(j))*RR(j,i))/v(2) + ((tbl(3)/R(j))*RR(j,i))/v(3) + (((R(j)-tbl(1)-tbl(2)-tbl(3))/R(j))*RR(j,i))/v(4);
        end
    end
end
TTA = TT;           %PURE TRAVEL TIMES
%NOISE INTERFERENCES
for i = 1:length(shot)
    for j = 1:length(R)
        a = TT(j,i)*randn(1)*w/100;
        if a < 0 && abs(a)>=TT(j,i)
            a = rand(1)*w/100;
        end
        TT(j,i) = TT(j,i) + a;
    end
end
%PLOTTING
figure(1)
L1x = [0 Offset Offset 0];L2x = L1x;L3x = L1x;L4x = L1x;
L1y = [0 0 bounds(1) bounds(1)];
L2y = [bounds(1) bounds(1) bounds(2) bounds(2)];
L3y = [bounds(2) bounds(2) bounds(3) bounds(3)];
L4y = [bounds(3) bounds(3) bounds(4) bounds(4)];

A = plot(L1x,L1y); fill(L1x,L1y,'c');hold on
B = plot(L2x,L2y); fill(L2x,L2y,'b');hold on
C = plot(L3x,L3y); fill(L3x,L3y,'g');hold on
D = plot(L4x,L4y); fill(L4x,L4y,'r');hold on

for p = 1:length(shot)
    plot(shot(p),0,'O','MarkerSize',5,'MarkerFaceColor','k');hold on   
end
for p = 1:length(R)
    plot(0,R(p),'V','MarkerSize',8,'MarkerFaceColor','y');hold on     
end

for i=1:length(shot)
    for j =1:length(R)
        plot([shot(i),0],[0,R(j)],'-k');hold on     %PLOTTING RAY TRACING
    end
end
title('FORWARD MODELLING DOWNHOLE')
ylabel('DEPTH(m)');xlabel('OFFSET(m)')
grid on;set(gca,'ydir','reverse')

%INVERSION ART 
S_guest = 0.0025;S_mod = ones(size(R))*S_guest;T=zeros(size(R));
n=10000;            %n = iterasi
Error = 1;toleransi=10^-7;i=0;pos_check = 1;
while i<n 
    for p=1:length(shot)
        for k = 1:length(R)
            if k==1
                b = R(k);a =(RR(k,p)*R(k)/b)^2;
                T(k) = S_mod(k)*RR(k,p);
                dt = TT(k,p)-T(k);
                dS1= dt*(RR(k,p)*(R(k)/b))/a;S_new1 =S_mod(k)+dS1;
                S_mod(k) = S_new1;
            elseif k==2
                b = R(k);a = ((RR(k,p)/b)*R(k-1))^2+((RR(k,p)/b)*(R(k)-R(k-1)))^2;
                T(k)    = (R(k-1)/b)*RR(k,p)*S_mod(k-1)+((R(k)-R(k-1))/b)*RR(k,p)*S_mod(k);
                dt      = TT(k,p)-T(k);
                dS1     = dt*((R(k-1)/b)*RR(k,p))/a;
                S_new1  = S_mod(k-1)+dS1;S_mod(k-1) = S_new1;
                dS2     = dt*(((R(k)-R(k-1))/b)*RR(k,p))/a;
                S_new2  = S_mod(k)+dS2;S_mod(k) = S_new2;
            elseif k==3
                b   = R(k);a =((RR(k,p)/b)*(R(k-2)))^2+((RR(k,p)/b)*(R(k-1)-R(k-2)))^2+((RR(k,p)/b)*(R(k)-R(k-1)))^2;
                T(k)=(RR(k,p)*R(k-2)/b)*S_mod(k-2)+(RR(k,p)*(R(k-1)-R(k-2))/b)*S_mod(k-1)+(RR(k,p)*(R(k)-R(k-1))/b)*S_mod(k);
                dt  = TT(k,p)-T(k);
                dS1 = dt*((R(k-2)/b)*RR(k,p))/a;S_new1=S_mod(k-2)+dS1;S_mod(k-2)=S_new1;
                dS2 = dt*(((R(k-1)-R(k-2))/b)*RR(k,p))/a;S_new2=S_mod(k-1)+dS2;S_mod(k-1)=S_new2;
                dS3 = dt*(((R(k)-R(k-1))/b)*RR(k,p))/a;S_new3=S_mod(k)+dS3;S_mod(k)=S_new3;
            elseif k==4
                b    = R(k);a=((RR(k,p)/b)*R(k-3))^2+((RR(k,p)/b)*(R(k-2)-R(k-3)))^2+((RR(k,p)/b)*(R(k-1)-R(k-2)))^2+((RR(k,p)/b)*(R(k)-R(k-1)))^2;
                T(k) =((RR(k,p)/b)*R(k-3))*S_mod(k-3)+((RR(k,p)/b)*(R(k-2)-R(k-3)))*S_mod(k-2)+((RR(k,p)/b)*(R(k-1)-R(k-2)))*S_mod(k-1)+((RR(k,p)/b)*(R(k)-R(k-1)))*S_mod(k);
                dt   = TT(k,p)-T(k);
                dS1  = dt*(R(k-3)/b)*RR(k,p)/a;S_new1=S_mod(k-3)+dS1;S_mod(k-3)=S_new1;
                dS2  = dt*((R(k-2)-R(k-3))/b)*RR(k,p)/a;S_new2=S_mod(k-2)+dS2;S_mod(k-2)=S_new2;
                dS3  = dt*((R(k-1)-R(k-2))/b)*RR(k,p)/a;S_new3=S_mod(k-1)+dS3;S_mod(k-1)=S_new3;
                dS4  = dt*((R(k)-R(k-1))/b)*RR(k,p)/a;S_new4=S_mod(k)+dS4;S_mod(k)=S_new4;
            elseif k==5
                b = R(k);a=((RR(k,p)/b)*R(k-4))^2+((RR(k,p)/b)*(R(k-3)-R(k-4)))^2+((RR(k,p)/b)*(R(k-2)-R(k-3)))^2+((RR(k,p)/b)*(R(k-1)-R(k-2)))^2+((RR(k,p)/b)*(R(k)-R(k-1)))^2;
                T(k) =(RR(k,p)/b)*(R(k-4))*S_mod(k-4)+(RR(k,p)/b)*(R(k-3)-R(k-4))*S_mod(k-3)+(RR(k,p)/b)*(R(k-2)-R(k-3))*S_mod(k-2)+(RR(k,p)/b)*(R(k-1)-R(k-2))*S_mod(k-1)+(RR(k,p)/b)*(R(k)-R(k-1))*S_mod(k);
                dt = TT(k,p)-T(k);
                dS1 = dt*((R(k-4)/b)*RR(k,p))/a;S_new1=S_mod(k-4)+dS1;S_mod(k-4)=S_new1;
                dS2 = dt*(((R(k-3)-R(k-4))/b)*RR(k,p))/a;S_new2=S_mod(k-3)+dS2;S_mod(k-3)=S_new2;
                dS3 = dt*(((R(k-2)-R(k-3))/b)*RR(k,p))/a;S_new3=S_mod(k-2)+dS3;S_mod(k-2)=S_new3;
                dS4 = dt*(((R(k-1)-R(k-2))/b)*RR(k,p))/a;S_new4=S_mod(k-1)+dS4;S_mod(k-1)=S_new4;
                dS5 = dt*((R(k)-R(k-1))/b)*RR(k,p)/a;S_new5=S_mod(k)+dS5;S_mod(k)=S_new5;
            else
                b =R(k);a=((RR(k,p)/b)*R(k-5))^2+((RR(k,p)/b)*(R(k-4)-R(k-5)))^2+((RR(k,p)/b)*(R(k-3)-R(k-4)))^2+((RR(k,p)/b)*(R(k-2)-R(k-3)))^2+((RR(k,p)/b)*(R(k-1)-R(k-2)))^2+((RR(k,p)/b)*(R(k)-R(k-1)))^2;
                T(k) =(RR(k,p)/b)*R(k-5)*S_mod(k-5)+(RR(k,p)/b)*(R(k-4)-R(k-5))*S_mod(k-4)+(RR(k,p)/b)*(R(k-3)-R(k-4))*S_mod(k-3)+(RR(k,p)/b)*(R(k-2)-R(k-3))*S_mod(k-2)+(RR(k,p)/b)*(R(k-1)-R(k-2))*S_mod(k-1)+(RR(k,p)/b)*(R(k)-R(k-1))*S_mod(k);
                dt = TT(k,p)-T(k);
                dS1 = dt*(R(k-5)/b)*RR(k,p)/a;S_new1=S_mod(k-5)+dS1;S_mod(k-5)=S_new1;
                dS2 = dt*((R(k-4)-R(k-5))/b)*RR(k,p)/a;S_new2=S_mod(k-4)+dS2;S_mod(k-4)=S_new2;
                dS3 = dt*((R(k-3)-R(k-4))/b)*RR(k,p)/a;S_new3=S_mod(k-3)+dS3;S_mod(k-3)=S_new3;
                dS4 = dt*((R(k-2)-R(k-3))/b)*RR(k,p)/a;S_new4=S_mod(k-2)+dS4;S_mod(k-2)=S_new4;
                dS5 = dt*((R(k-1)-R(k-2))/b)*RR(k,p)/a;S_new5=S_mod(k-1)+dS5;S_mod(k-1)=S_new5;
                dS6 = dt*((R(k)-R(k-1))/b)*RR(k,p)/a;S_new6=S_mod(k)+dS6;S_mod(k)=S_new6;
            end
        end
    end
    dt = TT(:,length(shot))-T';
    Error=(mean(sum((dt).^2)))^0.5;
    for l=1:length(S_mod)
        if S_mod(l)<0
            pos_check = 0;
        end
    end
    i = 1+i;
    if Error<=toleransi && pos_check == 1
        break
    end
end 

%PLOTTING INVERSI ART
figure(2)
S = S_mod;V = 1./S;V=V';V_mod=repmat(V,1,n_shot);
fprintf('ITERASI ART = %d\n',i)
fprintf('Kecepatan Inversi (1D)--->V\nWAKTU MODE INVERSI(1D)---T\nDELTA WAKTU INVERSI-FORWARD--->dt\n')
fprintf('Model Kecepatan Inversi(grid)--->V_mod\nTingkat Kesalahan--->Error\n')
V_add =V_mod(6,:);V_add_mod =vertcat(V_mod,V_add); 
Y_grid = [0 R];
[X,Y]= meshgrid(linspace(0,Offset,n_shot),Y_grid);
pcolor(X,Y,V_add_mod);hold on
colorbar;
set(gca,'Ydir','reverse');
title('Model Inversi Kecepatan ART');
xlabel('Offset Shot(m)');
ylabel('Depth Receiver(m)');
for p = 1:length(shot)
    plot(shot(p),0,'o','MarkerSize',5,'MarkerFaceColor','k');hold on
end
for p = 1:length(R)
    plot(0,R(p),'V','MarkerSize',8,'MarkerFaceColor','k');hold on
end
hold on

%SIRT INVERSION
S_Guest = 500;
Z = zeros(length(R),length(shot));Toleransi = 10^-7;ERROR = 1;
Z_mod = ones(size(R)).*S_Guest;j=0;iter=50000;Tt = zeros(size(T));
while j<iter 
    for p=1:length(shot)
        for k = 1:length(R)
            if k==1
                b = R(k);a =(RR(k,p)*R(k)/b)^2;
                Tt(k) = Z_mod(k)*RR(k,p);
                dt = TT(k,p)-Tt(k);
                dS1a= dt*(RR(k,p)*(R(k)/b))/a;
            elseif k==2
                b = R(k);a = ((RR(k,p)/b)*R(k-1))^2+((RR(k,p)/b)*(R(k)-R(k-1)))^2;
                Tt(k)    = (R(k-1)/b)*RR(k,p)*Z_mod(k-1)+((R(k)-R(k-1))/b)*RR(k,p)*Z_mod(k);
                dt      = TT(k,p)-Tt(k);
                dS1b     = dt*((R(k-1)/b)*RR(k,p))/a;
                dS2a     = dt*(((R(k)-R(k-1))/b)*RR(k,p))/a;
            elseif k==3
                b   = R(k);a =((RR(k,p)/b)*(R(k-2)))^2+((RR(k,p)/b)*(R(k-1)-R(k-2)))^2+((RR(k,p)/b)*(R(k)-R(k-1)))^2;
                Tt(k)=(RR(k,p)*R(k-2)/b)*Z_mod(k-2)+(RR(k,p)*(R(k-1)-R(k-2))/b)*Z_mod(k-1)+(RR(k,p)*(R(k)-R(k-1))/b)*Z_mod(k);
                dt  = TT(k,p)-Tt(k);
                dS1c = dt*((R(k-2)/b)*RR(k,p))/a;
                dS2b = dt*(((R(k-1)-R(k-2))/b)*RR(k,p))/a;
                dS3a = dt*(((R(k)-R(k-1))/b)*RR(k,p))/a;
            elseif k==4
                b    = R(k);a=((RR(k,p)/b)*R(k-3))^2+((RR(k,p)/b)*(R(k-2)-R(k-3)))^2+((RR(k,p)/b)*(R(k-1)-R(k-2)))^2+((RR(k,p)/b)*(R(k)-R(k-1)))^2;
                Tt(k) =((RR(k,p)/b)*R(k-3))*Z_mod(k-3)+((RR(k,p)/b)*(R(k-2)-R(k-3)))*Z_mod(k-2)+((RR(k,p)/b)*(R(k-1)-R(k-2)))*Z_mod(k-1)+((RR(k,p)/b)*(R(k)-R(k-1)))*Z_mod(k);
                dt   = TT(k,p)-Tt(k);
                dS1d  = dt*(R(k-3)/b)*RR(k,p)/a;
                dS2c  = dt*((R(k-2)-R(k-3))/b)*RR(k,p)/a;
                dS3b  = dt*((R(k-1)-R(k-2))/b)*RR(k,p)/a;
                dS4a  = dt*((R(k)-R(k-1))/b)*RR(k,p)/a;
            elseif k==5
                b = R(k);a=((RR(k,p)/b)*R(k-4))^2+((RR(k,p)/b)*(R(k-3)-R(k-4)))^2+((RR(k,p)/b)*(R(k-2)-R(k-3)))^2+((RR(k,p)/b)*(R(k-1)-R(k-2)))^2+((RR(k,p)/b)*(R(k)-R(k-1)))^2;
                Tt(k) =(RR(k,p)/b)*(R(k-4))*Z_mod(k-4)+(RR(k,p)/b)*(R(k-3)-R(k-4))*Z_mod(k-3)+(RR(k,p)/b)*(R(k-2)-R(k-3))*Z_mod(k-2)+(RR(k,p)/b)*(R(k-1)-R(k-2))*Z_mod(k-1)+(RR(k,p)/b)*(R(k)-R(k-1))*Z_mod(k);
                dt = TT(k,p)-Tt(k);
                dS1e = dt*((R(k-4)/b)*RR(k,p))/a;
                dS2d = dt*(((R(k-3)-R(k-4))/b)*RR(k,p))/a;
                dS3c = dt*(((R(k-2)-R(k-3))/b)*RR(k,p))/a;
                dS4b = dt*(((R(k-1)-R(k-2))/b)*RR(k,p))/a;
                dS5a = dt*((R(k)-R(k-1))/b)*RR(k,p)/a;
            else
                b =R(k);a=((RR(k,p)/b)*R(k-5))^2+((RR(k,p)/b)*(R(k-4)-R(k-5)))^2+((RR(k,p)/b)*(R(k-3)-R(k-4)))^2+((RR(k,p)/b)*(R(k-2)-R(k-3)))^2+((RR(k,p)/b)*(R(k-1)-R(k-2)))^2+((RR(k,p)/b)*(R(k)-R(k-1)))^2;
                Tt(k) =(RR(k,p)/b)*R(k-5)*Z_mod(k-5)+(RR(k,p)/b)*(R(k-4)-R(k-5))*Z_mod(k-4)+(RR(k,p)/b)*(R(k-3)-R(k-4))*Z_mod(k-3)+(RR(k,p)/b)*(R(k-2)-R(k-3))*Z_mod(k-2)+(RR(k,p)/b)*(R(k-1)-R(k-2))*Z_mod(k-1)+(RR(k,p)/b)*(R(k)-R(k-1))*Z_mod(k);
                dt = TT(k,p)-Tt(k);
                dS1f = dt*(R(k-5)/b)*RR(k,p)/a;
                dS2e = dt*((R(k-4)-R(k-5))/b)*RR(k,p)/a;
                dS3d = dt*((R(k-3)-R(k-4))/b)*RR(k,p)/a;
                dS4c = dt*((R(k-2)-R(k-3))/b)*RR(k,p)/a;
                dS5b = dt*((R(k-1)-R(k-2))/b)*RR(k,p)/a;
                dS6a = dt*((R(k)-R(k-1))/b)*RR(k,p)/a;
            end
        end
        dS1  =dS1a+dS1b+dS1c+dS1d+dS1e+dS1f;Z(1,p)= dS1;
        dS2  =dS2a+dS2b+dS2c+dS2d+dS2e;Z(2,p)=dS2;
        dS3  =dS3a+dS3b+dS3c+dS3d;Z(3,p)=dS3;
        dS4  =dS4a+dS4b+dS4c;Z(4,p)=dS4;
        dS5  =dS5a+dS5b;Z(5,p)=dS5;
        dS6  =dS6a;Z(6,p)=dS6;
    end
    Z_mod(1) = Z_mod(1)+mean(Z(1,:))/6;Z_mod(2) = Z_mod(2)+mean(Z(2,:))/5;Z_mod(3) = Z_mod(3)+mean(Z(3,:))/4;
    Z_mod(4) = Z_mod(4)+mean(Z(4,:))/3;Z_mod(5) = Z_mod(5)+mean(Z(5,:))/2;Z_mod(6) = Z_mod(6)+mean(Z(6,:));
    dt = TT(:,length(shot))-Tt';
    ERROR=(mean(sum((dt).^2)))^0.5;
    for l=1:length(Z_mod)
        if Z_mod(l)<0
            pos_check = 0;
        end
    end
    j = 1+j;
    if ERROR<=Toleransi && pos_check == 1
        break
    end
end

v = 1./Z_mod;

%PLOTTING INVERSI SIRT
figure(3)
S_SIRT = repmat(Z_mod',1,n_shot);v=v';v_mod=repmat(v,1,n_shot);
fprintf('ITERASI SIRT = %d\n',j)
fprintf('Kecepatan Inversi (1D)--->v\nWAKTU MODE INVERSI(1D)---Tt\nDELTA WAKTU INVERSI-FORWARD--->dt\n')
fprintf('Model Kecepatan Inversi(grid)--->v_mod\nTingkat Kesalahan--->ERROR\n')
v_add =v_mod(6,:);v_add_mod =vertcat(v_mod,v_add); 
Y_grid = [0 R];
[X,Y]= meshgrid(linspace(0,Offset,n_shot),Y_grid);
pcolor(X,Y,v_add_mod);hold on
colorbar;
set(gca,'Ydir','reverse');
title('Model Inversi Kecepatan SIRT');
xlabel('Offset Shot(m)');
ylabel('Depth Receiver(m)');
for p = 1:length(shot)
    plot(shot(p),0,'o','MarkerSize',5,'MarkerFaceColor','k');hold on
end
for p = 1:length(R)
    plot(0,R(p),'V','MarkerSize',8,'MarkerFaceColor','k');hold on
end
hold on


