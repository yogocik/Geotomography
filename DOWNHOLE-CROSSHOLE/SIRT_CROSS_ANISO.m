%UAS GEOTOMOGRAPHY NO 1 (MODEL 3 LAPIS 7 REC-SHOT)
%KELOMPOK 1 
%INVERSI MEDIUM ANISOTROPI METODE SIRT MODEL CROSSHOLE 
%VARIASI PADA SIFAT ANISOTROPI LAPISAN 2&3 serta KETEBALAN PER LAPISAN
clc;clear;close all
Depth = 30;Offset = Depth;                  %Dimensi Area
vlap = [200,300,400];
slap = 1./vlap;                             %Slowness
tbl = [20,5,5];                             %Tebal Lapisan
bounds = [tbl(1) sum(tbl(1:2)) sum(tbl)];   %Batas Lapisan
dlap = [0,2,1];eps=2;                       %Parameter Anisotropi   
int.rec = 4.5; int.shot = int.rec;          %Offset Shot-Receiver
rec = 0:int.rec:Depth;shot = rec;
w = 0;                                      %Random Noise Ratio
iter=1000;                                  %ITERASI INVERSI SIRT
fprintf('LAPISAN 1 ISOTROPI (PERMANENT-DEFAULT)\n')
if dlap(2)~=0
    fprintf('LAPISAN 2 ANISOTROPI\n')
else
    fprintf('LAPISAN 2 ISOTROPI\n')
end
if dlap(3)~=0
    fprintf('LAPISAN 3 ANISOTROPI\n')
else
    fprintf('LAPISAN 3 ISOTROPI\n')
end
if bounds(3)~=Depth
    error('Total tebal lapisan tidak sama dengan Depth (input bounds)')
end
for i=1:length(tbl)
    if mod(tbl(i),5)~=0
        error('Tebal lapisan harus minimal 5 dan kelipatan 5 (input tbl)')
    end
end
if dlap(1)~=0
    error('Lapisan pertama harus isotropi (dlap(1) = 0)')
end
if rec(1)>= bounds(1)
    error('Harus ada shot-receiver di lapisan pertama (tidak boleh pada batas lapisan 1)')
end
count = 0;
for i = 1:length(dlap)
    if dlap(i)==0
        count = count + 1;
    end
end
if count == 3
    error('HARUS ADA MINIMAL SATU LAPISAN ANISOTROPI PADA LAPISAN 2 atau 3 (input nilai dlap)')
end
RP = zeros(length(rec),length(shot));teta = zeros(size(RP));
for i=1:length(shot)
    for j=1:length(rec)
        RP(j,i) = sqrt(Offset^2+(rec(j)-shot(i))^2);
        teta(j,i)= 90-acosd(Offset/RP(j,i));
    end
end
figure(1)
L1x = [0 Offset Offset 0];L2x = L1x;L3x = L1x;L4x = L1x;L5x = L1x;
L1y = [0 0 bounds(1) bounds(1)];
L2y = [bounds(1) bounds(1) bounds(2) bounds(2)];
L3y = [bounds(2) bounds(2) bounds(3) bounds(3)];
A = plot(L1x,L1y); fill(L1x,L1y,'c');hold on
B = plot(L2x,L2y); fill(L2x,L2y,'b');hold on
C = plot(L3x,L3y); fill(L3x,L3y,'g');hold on

for i = 1:length(shot)
    for j = 1:length(rec)
        line([0 Offset],[shot(i) rec(j)],'Color','k','LineStyle','-','LineWidth',1);hold on
    end
    plot(0,shot(i),'v','MarkerSize',8,'MarkerFaceColor','r');hold on
    plot(Offset,rec(i),'o','MarkerSize',8,'MarkerFaceColor','r');hold on
end
title('CROSS-HOLE ANISOTROPY')
xlabel('Offset(m)');ylabel('Depth(m)')
set(gca,'Ydir','reverse')
TT = zeros(size(RP));
for i = 1:length(shot)
    for j = 1:length(rec)
        if shot(i)<bounds(1)
            if rec(j)<bounds(1)
                TT(j,i) = RP(j,i)*slap(1);
            elseif rec(j)<bounds(2) && rec(j)>bounds(1)
                if dlap(2)~=0
                    p1 = sind(teta(j,i))*slap(1);
                    y1 = [p1*vlap(2)*(eps-dlap(2)) 0 p1*vlap(2)*dlap(2) -1 p1*vlap(2)];
                    k1 = abs(asind(roots(y1)));
                    if length(k1) > 1
                        for jj = 1:length(k1)
                            k11 = k1(jj)-teta(j,i);
                        end
                        [val,index] = min(k11);k1 = k1(index);
                    end
                    TT(j,i) = ((bounds(1)-shot(i))/(rec(j)-shot(i)))*RP(j,i)*slap(1) + ((rec(j)-bounds(1))/(rec(j)-shot(i)))*RP(j,i)*slap(2)*(1+(-1*sind(k1)^2*cosd(k1)^2*dlap(2)));
                end
                TT(j,i) = ((bounds(1)-shot(i))/(rec(j)-shot(i)))*RP(j,i)*slap(1)+ ((rec(j)-bounds(1))/(rec(j)-shot(i)))*RP(j,i)*slap(2);
            elseif rec(j)<bounds(3) && rec(j)>bounds(2)
                if dlap(2)~=0 && dlap(3)~=0
                    p2 = sind(teta(j,i))*slap(2)*(1+(-1*sind(k1)^2*cosd(k1)^2*dlap(2)));
                    y2 = [p2*vlap(3)*(eps-dlap(3)) 0 p2*vlap(3)*dlap(3) -1 p2*vlap(3)];
                    k2 = abs(asind(roots(y2)));
                    if length(k2) > 1
                        for jj = 1:length(k1)
                            k22 = k2(jj)-teta(j,i);
                        end
                        [val,index] = min(k22);k2 = k2(index);
                    end
                    TT(j,i) = ((bounds(1)-shot(i))/(rec(j)-shot(i)))*RP(j,i)*slap(1) + ((bounds(2)-bounds(1))/(rec(j)-shot(i)))*RP(j,i)*slap(2)*(1+(-1*sind(k1)^2*cosd(k1)^2*dlap(2)))+...
                        ((rec(j)-bounds(2))/(rec(j)-shot(i)))*RP(j,i)*slap(3)*(1+(-1*sind(k2)^2*cosd(k2)^2*dlap(3)));
                elseif dlap(2)==0 && dlap(3)~=0
                    p2 = sind(teta(j,i))*slap(2);
                    y2 = [p2*vlap(3)*(eps-dlap(3)) 0 p2*vlap(3)*dlap(3) -1 p2*vlap(3)];
                    k2 = abs(asind(roots(y2)));
                    if length(k2) > 1
                        for jj = 1:length(k2)
                            k22 = k2(jj)-teta(j,i);
                        end
                        [val,index] = min(k22);k2 = k2(index);
                    end
                    TT(j,i) = ((bounds(1)-shot(i))/(rec(j)-shot(i)))*RP(j,i)*slap(1) + ((bounds(2)-bounds(1))/(rec(j)-shot(i)))*RP(j,i)*slap(2)+...
                        ((rec(j)-bounds(2))/(rec(j)-shot(i)))*RP(j,i)*slap(3)*(1+-1*sind(k2)^2*cosd(k2)^2*dlap(3));
                elseif dlap(2) ~= 0 && dlap(3) == 0
                    p1 = sind(teta(j,i))*slap(1);
                    y1 = [p1*vlap(2)*(eps-dlap(2)) 0 p1*vlap(2)*dlap(2) -1 p1*vlap(2)];
                    k1 = abs(asind(roots(y1)));
                    if length(k1) > 1
                        for jj = 1:length(k1)
                            k11 = k1(jj)-teta(j,i);
                        end
                        [val,index] = min(k11);k1 = k1(index);
                    end
                    TT(j,i) = ((bounds(1)-shot(i))/(rec(j)-shot(i)))*RP(j,i)*slap(1) + ((rec(j)-bounds(1))/(rec(j)-shot(i)))*RP(j,i)*slap(2)*(1+(-1*sind(k1)^2*cosd(k1)^2*dlap(2)))+...
                        ((rec(j)-bounds(2))/(rec(j)-shot(i)))*RP(j,i)*slap(3);
                end
            end
        elseif shot(i)<bounds(2)
            if rec(j)<bounds(1)
                TT(j,i) = ((shot(i)-bounds(1))/shot(i))*RP(j,i)*slap(2)+(bounds(1)/shot(i))*RP(j,i)*slap(1);
            elseif rec(j)<bounds(2)
                TT(j,i) = RP(j,i)*slap(2);
            elseif rec(j)<bounds(3)
                if dlap(3)~=0
                    p3 = sind(teta(j,i))*slap(2);
                    y3 = [p3*vlap(3)*(eps-dlap(3)) 0 p2*vlap(3)*dlap(3) -1 p2*vlap(3)];
                    k3 = abs(asind(roots(y3)));
                    if length(k3) > 1
                        for jj = 1:length(k3)
                            k33 = k3(jj)-teta(j,i);
                        end
                        [val,index] = min(k33);k3 = k3(index);
                    end
                    TT(j,i) = (bounds(2)-shot(i))/(rec(j)-shot(i))*RP(j,i)*slap(2)+((rec(j)-bounds(2))/(rec(j)-shot(i)))*RP(j,i)*slap(3)*(1+(-1)*sind(k3)^2*cosd(k3)^2*dlap(3));
                end
                TT(j,i) = (bounds(2)-shot(i))/(rec(j)-shot(i))*RP(j,i)*slap(2)+((rec(j)-bounds(2))/(rec(j)-shot(i)))*RP(j,i)*slap(3);
            end
        elseif shot(i)<bounds(3)
            if rec(j)<bounds(1)
                if dlap(2) ~= 0
                    p32 = sind(teta(j,i))*slap(3);
                    y32 = [p32*vlap(2)*(eps-dlap(2)) 0 p32*vlap(2)*dlap(2) -1 p32*vlap(2)];
                    k32 = abs(asind(roots(y32)));
                    if length(k32) > 1
                        for jj = 1:length(k32)
                            k332 = k32(jj)-teta(j,i);
                        end
                        [val,index] = min(k332);k32 = k32(index);
                    end
                    TT(j,i)=(shot(i)-bounds(2))/(shot(i)-rec(j))*RP(j,i)*slap(3)+((bounds(2)-bounds(1))/(shot(i)-rec(j)))*RP(j,i)*slap(2)*(1+(-1)*sind(k32)^2*cosd(k32)^2*dlap(2))+...
                        ((bounds(1)-rec(j))/(shot(i)-rec(j)))*RP(j,i)*slap(1);
                end
                TT(j,i) = (shot(i)-bounds(2))/(shot(i)-rec(j))*RP(j,i)*slap(3)+((bounds(2)-bounds(1))/(shot(i)-rec(j)))*RP(j,i)*slap(2)+...
                        ((bounds(1)-rec(j))/(shot(i)-rec(j)))*RP(j,i)*slap(1);
            elseif rec(j)<bounds(2)
                if dlap(2)~=0
                    TT(j,i)=((shot(i)-bounds(2))/(shot(i)-rec(j)))*RP(j,i)*slap(3)+((bounds(2)-rec(j))/(shot(i)-rec(j)))*RP(j,i)*slap(2)*(1+(-1)*sind(k32)^2*cosd(k32)^2*dlap(2));
                end
                TT(j,i) = ((shot(i)-bounds(2))/(shot(i)-rec(j)))*RP(j,i)*slap(3)+((bounds(2)-rec(j))/(shot(i)-rec(j)))*RP(j,i)*slap(2);
            elseif rec(j)<bounds(3)
                TT(j,i) = RP(j,i)*slap(3);
            end
        end
    end
end
TTA = TT;           %PURE TRAVEL TIMES
%NOISE INTERFERENCES
for i = 1:length(shot)
    for j = 1:length(rec)
        a = TT(j,i)*randn(1)*w/100;
        if a < 0 && abs(a)>=TT(j,i)
            a = rand(1)*w/100;
        end
        TT(j,i) = TT(j,i) + a;
    end
end
fprintf('TRAVEL TIME NO-NOISE : TTA\n')
fprintf('RANDOM-NOISE RATIO : %d persen\n',w)
fprintf('TRAVEL TIME WITH-NOISE : TT\n')

%SIRT INVERSION
S_Guest = 300; AA = 1;
Z = zeros(length(rec),length(shot));Toleransi = 10^-2;ERROR = 1;
Z_mod = ones(length(rec),1).*S_Guest;p=0;Tt = zeros(size(Z));D = zeros(size(Tt));
A_mod = ones(size(Z_mod)).*AA;g = zeros(size(A_mod));PERTU_S = zeros(size(Z_mod));PERTU_A=PERTU_S;
dS = zeros(size(Z_mod));dA=zeros(size(dS));TOLERANSI_UP = 10^5;
%g = perubahan slowness aniso akibat teta,
%A_mod adalah Delta Aniso Tebakan, D = Selisih Travel Time Inversi-Forward
while p<iter
    for i = 1:length(shot)
        for j = 1:length(rec)
            if shot(i) == shot(1)
                if rec(j) == shot(i)
                    continue
                end
                if rec(j) == shot(2)
                    TETA(1,1) = 0;g(1,1) = 0; 
                    Tt(j,i) = RP(j,i)*(Z_mod(1,1)+g(1,1)*A_mod(1,1));
                    dt   = TT(j,i)-Tt(j,i);D(j,i) = dt;
                    dS1b = dt*RP(j,i)/((1+g(1,1)^2)*RP(j,i)^2);
                    dA1b = dt*RP(j,i)*g(1,1)/(1+g(1,1)^2)*RP(j,i)^2;
                elseif rec(j) == shot(3)
                    p0 = sind(teta(j,i))*(Z_mod(1,1));
                    y0 = [p0*(1/Z_mod(2,1))*(eps-A_mod(2,1)) 0 p0*(1/Z_mod(2,1))*A_mod(2,1) -1 p0*(1/Z_mod(2,1))];
                    k0 = abs(asind(roots(y0)));
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(2,1)  = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(2)-rec(1))/rec(3))*RP(j,i)*(Z_mod(1,1)+g(1,1)*A_mod(1,1))+...
                        ((rec(3)-rec(2))/rec(3))*RP(j,i)*(Z_mod(2,1)+ g(2,1)*A_mod(2,1));
                    dt      = TT(j,i) - Tt(j,i); D(j,i) = dt;
                    dd = ((rec(2)-rec(1))^2+(rec(3)-rec(2))^2)*(RP(j,i)/rec(3))^2;
                    dS1c = dt*((rec(2)-rec(1))/rec(3))*RP(j,i)/((1+g(1,1)^2)*dd);
                    dA1c = dt*((rec(2)-rec(1))/rec(3))*RP(j,i)*g(1,1)/((1+g(1,1)^2)*dd);
                    dS2a = dt*((rec(3)-rec(2))/rec(3))*RP(j,i)/((1+g(2,1)^2)*dd);                     %grid 2 Vertikal (antara shot 2 dan 3)
                    dA2a = dt*((rec(3)-rec(2))/rec(3))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                elseif rec(j) == shot(4)
                    p0 = sind(teta(j,i))*(Z_mod(1,1)+g(1,1)*A_mod(1,1));
                    y0 = [p0*(1/Z_mod(2,1))*(eps-A_mod(2,1)) 0 p0*(1/Z_mod(2,1))*A_mod(2,1) -1 p0*(1/Z_mod(2,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(2,1)  = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(2,1)+g(2,1)*A_mod(2,1));
                    y0 = [p0*(1/Z_mod(3,1))*(eps-A_mod(3,1)) 0 p0*(1/Z_mod(3,1))*A_mod(3,1) -1 p0*(1/Z_mod(3,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1)  = -1*sind(teta(j,i))^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(2)-rec(1))/rec(4))*RP(j,i)*(Z_mod(1,1)+g(1,1)*A_mod(1,1))+...
                        ((rec(3)-rec(2))/rec(4))*RP(j,i)*(Z_mod(2,1)+g(2,1)*A_mod(2,1))+...
                        ((rec(4)-rec(3))/rec(4))*RP(j,i)*(Z_mod(3,1)+g(3,1)*A_mod(3,1));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd = ((rec(2)-rec(1))^2+(rec(3)-rec(2))^2+(rec(4)-rec(3))^2)*(RP(j,i)/rec(4))^2;
                    dS1d = dt*((rec(2)-rec(1))/rec(4))*RP(j,i)/((1+g(1,1)^2)*dd);
                    dA1d = dt*((rec(2)-rec(1))/rec(4))*RP(j,i)*g(1,1)/((1+g(1,1)^2)*dd);
                    dS2b = dt*((rec(3)-rec(2))/rec(4))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA2b = dt*((rec(3)-rec(2))/rec(4))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS3a = dt*((rec(4)-rec(3))/rec(4))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA3a = dt*((rec(4)-rec(3))/rec(4))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                elseif rec(j) == shot(5)
                    p0 = sind(teta(j,i))*(Z_mod(1,1)+g(1,1)*A_mod(1,1));
                    y0 = [p0*(1/Z_mod(2,1))*(eps-A_mod(2,1)) 0 p0*(1/Z_mod(2,1))*A_mod(2,1) -1 p0*(1/Z_mod(2,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(2,1)  = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(2,1)+g(2,1)*A_mod(2,1));
                    y0 = [p0*(1/Z_mod(3,1))*(eps-A_mod(3,1)) 0 p0*(1/Z_mod(3,1))*A_mod(3,1) -1 p0*(1/Z_mod(3,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1)  = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(3,1)+g(3,1)*A_mod(3,1));
                    y0 = [p0*(1/Z_mod(4,1))*(eps-A_mod(4,1)) 0 p0*(1/Z_mod(4,1))*A_mod(4,1) -1 p0*(1/Z_mod(4,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(2)-rec(1))/rec(5))*RP(j,i)*(Z_mod(1,1)+g(1,1)*A_mod(1,1))+...
                        ((rec(3)-rec(2))/rec(5))*RP(j,i)*(Z_mod(2,1)+g(2,1)*A_mod(2,1))+...
                        ((rec(4)-rec(3))/rec(5))*RP(j,i)*(Z_mod(3,1)+g(3,1)*A_mod(3,1))+...
                        ((rec(5)-rec(4))/rec(5))*RP(j,i)*(Z_mod(4,1)+g(4,1)*A_mod(4,1));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd = ((rec(2)-rec(1))^2+(rec(3)-rec(2))^2+(rec(4)-rec(3))^2+(rec(5)-rec(4))^2)*(RP(j,i)/rec(5))^2;
                    dS1e = dt*((rec(2)-rec(1))/rec(5))*RP(j,i)/((1+g(1,1)^2)*RP(j,i)^2);
                    dA1e = dt*((rec(2)-rec(1))/rec(5))*RP(j,i)*g(1,1)/((1+g(1,1)^2)*dd);
                    dS2c = dt*((rec(3)-rec(2))/rec(5))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA2c = dt*((rec(3)-rec(2))/rec(5))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS3b = dt*((rec(4)-rec(3))/rec(5))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA3b = dt*((rec(4)-rec(3))/rec(5))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS4a = dt*((rec(5)-rec(4))/rec(5))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA4a = dt*((rec(5)-rec(4))/rec(5))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                elseif rec(j) == shot(6)
                    p0 = sind(teta(j,i))*(Z_mod(1,1)+g(1,1)*A_mod(1,1));
                    y0 = [p0*(1/Z_mod(2,1))*(eps-A_mod(2,1)) 0 p0*(1/Z_mod(2,1))*A_mod(2,1) -1 p0*(1/Z_mod(2,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(2,1)  = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(2,1)+g(2,1)*A_mod(2,1));
                    y0 = [p0*(1/Z_mod(3,1))*(eps-A_mod(3,1)) 0 p0*(1/Z_mod(3,1))*A_mod(3,1) -1 p0*(1/Z_mod(3,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1)  = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(3,1)+g(3,1)*A_mod(3,1));
                    y0 = [p0*(1/Z_mod(4,1))*(eps-A_mod(4,1)) 0 p0*(1/Z_mod(4,1))*A_mod(4,1) -1 p0*(1/Z_mod(4,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(4,1)+g(4,1)*A_mod(4,1));
                    y0 = [p0*(1/Z_mod(5,1))*(eps-A_mod(5,1)) 0 p0*(1/Z_mod(5,1))*A_mod(5,1) -1 p0*(1/Z_mod(5,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(5,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(2)-rec(1))/rec(6))*RP(j,i)*(Z_mod(1,1)+g(1,1)*A_mod(1,1))+...
                        ((rec(3)-rec(2))/rec(6))*RP(j,i)*(Z_mod(2,1)+g(2,1)*A_mod(2,1))+...
                        ((rec(4)-rec(3))/rec(6))*RP(j,i)*(Z_mod(3,1)+g(3,1)*A_mod(3,1))+...
                        ((rec(5)-rec(4))/rec(6))*RP(j,i)*(Z_mod(4,1)+g(4,1)*A_mod(4,1))+...
                        ((rec(6)-rec(5))/rec(6))*RP(j,i)*(Z_mod(5,1)+g(5,1)*A_mod(5,1));
                    dt = TT(j,i)- Tt(j,i); D(j,i) = dt;
                    dd = ((rec(2)-rec(1))^2+(rec(3)-rec(2))^2+(rec(4)-rec(3))^2+(rec(5)-rec(4))^2+(rec(6)-rec(5))^2)*(RP(j,i)/rec(6))^2; 
                    dS1f = dt*((rec(2)-rec(1))/rec(6))*RP(j,i)/((1+g(1,1)^2)*dd);
                    dA1f = dt*((rec(2)-rec(1))/rec(6))*RP(j,i)*g(1,1)/((1+g(1,1)^2)*dd);
                    dS2d = dt*((rec(3)-rec(2))/rec(6))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA2d = dt*((rec(3)-rec(2))/rec(6))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS3c = dt*((rec(4)-rec(3))/rec(6))*RP(j,i)/((1+g(3,1)^2)*RP(j,i)^2);
                    dA3c = dt*((rec(4)-rec(3))/rec(6))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS4b = dt*((rec(5)-rec(4))/rec(6))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA4b = dt*((rec(5)-rec(4))/rec(6))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS5a = dt*((rec(6)-rec(5))/rec(6))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA5a = dt*((rec(6)-rec(5))/rec(6))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                elseif rec(j) == shot(7)
                    p0 = sind(teta(j,i))*(Z_mod(1,1)+g(1,1)*A_mod(1,1));
                    y0 = [p0*(1/Z_mod(2,1))*(eps-A_mod(2,1)) 0 p0*(1/Z_mod(2,1))*A_mod(2,1) -1 p0*(1/Z_mod(2,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(2,1)  = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(2,1)+g(2,1)*A_mod(2,1));
                    y0 = [p0*(1/Z_mod(3,1))*(eps-A_mod(3,1)) 0 p0*(1/Z_mod(3,1))*A_mod(3,1) -1 p0*(1/Z_mod(3,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1)  = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(3,1)+g(3,1)*A_mod(3,1));
                    y0 = [p0*(1/Z_mod(4,1))*(eps-A_mod(4,1)) 0 p0*(1/Z_mod(4,1))*A_mod(4,1) -1 p0*(1/Z_mod(4,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(4,1)+g(4,1)*A_mod(4,1));
                    y0 = [p0*(1/Z_mod(5,1))*(eps-A_mod(5,1)) 0 p0*(1/Z_mod(5,1))*A_mod(5,1) -1 p0*(1/Z_mod(5,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(5,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(5,1)+g(5,1)*A_mod(5,1));
                    y0 = [p0*(1/Z_mod(6,1))*(eps-A_mod(6,1)) 0 p0*(1/Z_mod(6,1))*A_mod(6,1) -1 p0*(1/Z_mod(6,1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(6,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(2)-rec(1))/rec(7))*RP(j,i)*(Z_mod(1,1)+g(1,1)*A_mod(1,1))+...
                        ((rec(3)-rec(2))/rec(7))*RP(j,i)*(Z_mod(2,1)+g(2,1)*A_mod(2,1))+...
                        ((rec(4)-rec(3))/rec(7))*RP(j,i)*(Z_mod(3,1)+g(3,1)*A_mod(3,1))+...
                        ((rec(5)-rec(4))/rec(7))*RP(j,i)*(Z_mod(4,1)+g(4,1)*A_mod(4,1))+...
                        ((rec(6)-rec(5))/rec(7))*RP(j,i)*(Z_mod(5,1)+g(5,1)*A_mod(5,1))+...
                        ((rec(7)-rec(6))/rec(7))*RP(j,i)*(Z_mod(6,1)+g(6,1)*A_mod(6,1));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd = ((rec(2)-rec(1))^2+(rec(3)-rec(2))^2+(rec(4)-rec(3))^2+(rec(5)-rec(4))^2+(rec(6)-rec(5))^2+(rec(7)-rec(6))^2)*(RP(j,i)/rec(7))^2;
                    dS1g = dt*((rec(2)-rec(1))/rec(7))*RP(j,i)/((1+g(1,1)^2)*dd);
                    dA1g = dt*((rec(2)-rec(1))/rec(7))*RP(j,i)*g(1,1)/((1+g(1,1)^2)*dd);
                    dS2e = dt*((rec(3)-rec(2))/rec(7))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA2e = dt*((rec(3)-rec(2))/rec(7))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS3d = dt*((rec(4)-rec(3))/rec(7))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA3d = dt*((rec(4)-rec(3))/rec(7))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS4c = dt*((rec(5)-rec(4))/rec(7))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA4c = dt*((rec(5)-rec(4))/rec(7))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS5b = dt*((rec(6)-rec(5))/rec(7))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA5b = dt*((rec(6)-rec(5))/rec(7))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                    dS6a = dt*((rec(7)-rec(6))/rec(7))*RP(j,i)/((1+g(6,1)^2)*dd);
                    dA6a = dt*((rec(7)-rec(6))/rec(7))*RP(j,i)*g(6,1)/((1+g(6,1)^2)*dd);
                    %RATA SIRT TEMBAKAN
                    dS1A = dS1b+dS1c+dS1d+dS1e+dS1f+dS1g;dS(1,i) = dS1A;
                    dA1A = dA1b+dA1c+dA1d+dA1e+dA1f+dA1g;dA(1,i) = dA1A;
                    dS2A = dS2a+dS2b+dS2c+dS2d+dS2e;dS(2,i) = dS2A;
                    dA2A = dA2a+dA2b+dA2c+dA2d+dA2e;dA(2,i) = dA2A;
                    dS3A = dS3a+dS3b+dS3c+dS3d;dS(3,i)=dS3A;
                    dA3A = dA3a+dA3b+dA3c+dA3d;dA(3,i)=dA3A;
                    dS4A = dS4a+dS4b+dS4c;dS(4,i)=dS4A;
                    dA4A = dA4a+dA4b+dA4c;dA(4,i)=dA4A;
                    dS5A = dS5a+dS5b;dS(5,i)=dS5A;
                    dA5A = dA5a+dA5b;dA(5,i)=dA5A;
                    dS6A = dS6a;dS(6,i)=dS6A;
                    dA6A = dA6a;dA(6,i)=dA6A;
                end
            elseif shot(i) == shot(2)
                if rec(j) == shot(i)
                    continue
                end
                if rec(j) == shot(1)
                    Tt(j,i) = RP(j,i)*(Z_mod(1)+g(1,1)*A_mod(1));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dS11a = dt*((rec(2)-rec(1))/rec(2))*RP(j,i)/((1+g(1,1)^2)*RP(j,i)^2);
                    dA11a = dt*((rec(2)-rec(1))/rec(2))*RP(j,i)*g(1,1)/((1+g(1,1)^2)*RP(j,i)^2);
                elseif rec(j) == shot(3)
                    g(2,1) = g(1,1);
                    TT(j,i) = RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2));
                    dt = TT(j,i)-Tt(j,i);D(j,i) = dt;
                    dS22a = dt*((rec(3)-rec(2))/rec(2))*RP(j,i)/((1+g(2,1)^2)*RP(j,i)^2);
                    dA22a = dt*((rec(3)-rec(2))/rec(2))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*RP(j,i)^2);
                elseif rec(j) == shot(4)
                    p0 = sind(teta(j,i))*(Z_mod(2)+g(2,1)*A_mod(2));
                    y0 = [p0*(1/Z_mod(3))*(eps-A_mod(3)) 0 p0*(1/Z_mod(3))*A_mod(3) -1 p0*(1/Z_mod(3))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(3)-rec(2))/(rec(4)-rec(2)))*RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2))+...
                        ((rec(4)-rec(3))/(rec(4)-rec(2)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd = ((rec(3)-rec(2))^2+(rec(4)-rec(3))^2)*(RP(j,i)/(rec(4)-rec(2)))^2;
                    dS22b = dt*((rec(3)-rec(2))/(rec(4)-rec(2)))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA22b = dt*((rec(3)-rec(2))/(rec(4)-rec(2)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS33a = dt*((rec(4)-rec(3))/(rec(4)-rec(2)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA33a = dt*((rec(4)-rec(3))/(rec(4)-rec(2)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                elseif rec(j) == shot(5)
                    p0 = sind(teta(j,i))*(Z_mod(2)+g(2,1)*A_mod(2));
                    y0 = [p0*(1/Z_mod(3))*(eps-A_mod(3)) 0 p0*(1/Z_mod(3))*A_mod(3) -1 p0*(1/Z_mod(3))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(3)-rec(2))/(rec(5)-rec(2)))*RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2))+...
                        ((rec(4)-rec(3))/(rec(5)-rec(2)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(5)-rec(2)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(3)-rec(2))^2+(rec(4)-rec(3))^2+(rec(5)-rec(4))^2)*(RP(j,i)/(rec(5)-rec(2)))^2;
                    dS22c = dt*((rec(3)-rec(2))/(rec(5)-rec(2)))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA22c = dt*((rec(3)-rec(2))/(rec(5)-rec(2)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS33b = dt*((rec(4)-rec(3))/(rec(5)-rec(2)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA33b = dt*((rec(4)-rec(3))/(rec(5)-rec(2)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS44a = dt*((rec(5)-rec(4))/(rec(5)-rec(2)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA44a = dt*((rec(5)-rec(4))/(rec(5)-rec(2)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                elseif rec(j) == shot(6)
                    p0 = sind(teta(j,i))*(Z_mod(2)+g(2,1)*A_mod(2));
                    y0 = [p0*(1/Z_mod(3))*(eps-A_mod(3)) 0 p0*(1/Z_mod(3))*A_mod(3) -1 p0*(1/Z_mod(3))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(5))*(eps-A_mod(5)) 0 p0*(1/Z_mod(5))*A_mod(5) -1 p0*(1/Z_mod(5))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(5,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(3)-rec(2))/(rec(6)-rec(2)))*RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2))+...
                        ((rec(4)-rec(3))/(rec(6)-rec(2)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(6)-rec(2)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(6)-rec(2)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(3)-rec(2))^2+(rec(4)-rec(3))^2+(rec(5)-rec(4))^2+(rec(6)-rec(5))^2)*(RP(j,i)/(rec(6)-rec(2)))^2;
                    dS22d = dt*((rec(3)-rec(2))/(rec(6)-rec(2)))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA22d = dt*((rec(3)-rec(2))/(rec(6)-rec(2)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS33c = dt*((rec(4)-rec(3))/(rec(6)-rec(2)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA33c = dt*((rec(4)-rec(3))/(rec(6)-rec(2)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS44b = dt*((rec(5)-rec(4))/(rec(6)-rec(2)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA44b = dt*((rec(5)-rec(4))/(rec(6)-rec(2)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS55a = dt*((rec(6)-rec(5))/(rec(6)-rec(2)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA55a = dt*((rec(6)-rec(5))/(rec(6)-rec(2)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                elseif rec(j) == shot(7)
                    p0 = sind(teta(j,i))*(Z_mod(2)+g(2,1)*A_mod(2));
                    y0 = [p0*(1/Z_mod(3))*(eps-A_mod(3)) 0 p0*(1/Z_mod(3))*A_mod(3) -1 p0*(1/Z_mod(3))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(5))*(eps-A_mod(5)) 0 p0*(1/Z_mod(5))*A_mod(5) -1 p0*(1/Z_mod(5))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(5,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(5)+g(5,1)*A_mod(5));
                    y0 = [p0*(1/Z_mod(6))*(eps-A_mod(6)) 0 p0*(1/Z_mod(6))*A_mod(6) -1 p0*(1/Z_mod(6))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(6,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(3)-rec(2))/(rec(7)-rec(2)))*RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2))+...
                        ((rec(4)-rec(3))/(rec(7)-rec(2)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(7)-rec(2)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(7)-rec(2)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(3)-rec(2))^2+(rec(4)-rec(3))^2+(rec(5)-rec(4))^2+(rec(6)-rec(5))^2+(rec(7)-rec(6))^2)*(RP(j,i)/(rec(7)-rec(2)))^2;
                    dS22e = dt*((rec(3)-rec(2))/(rec(7)-rec(2)))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA22e = dt*((rec(3)-rec(2))/(rec(7)-rec(2)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS33d = dt*((rec(4)-rec(3))/(rec(7)-rec(2)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA33d = dt*((rec(4)-rec(3))/(rec(7)-rec(2)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS44c = dt*((rec(5)-rec(4))/(rec(7)-rec(2)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA44c = dt*((rec(5)-rec(4))/(rec(7)-rec(2)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS55b = dt*((rec(6)-rec(5))/(rec(7)-rec(2)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA55b = dt*((rec(6)-rec(5))/(rec(7)-rec(2)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                    dS66a = dt*((rec(6)-rec(5))/(rec(7)-rec(2)))*RP(j,i)/((1+g(6,1)^2)*dd);
                    dA66a = dt*((rec(6)-rec(5))/(rec(7)-rec(2)))*RP(j,i)*g(6,1)/((1+g(6,1)^2)*dd);
                    %RATA SIRT TEMBAKAN
                    dS11 = dS11a;dS(1,i)=dS11;
                    dA11 = dA11a;dA(1,i)=dA11;
                    dS22 = dS22a+dS22b+dS22c+dS22d+dS22e;dS(2,i)=dS22;
                    dA22 = dA22a+dA22b+dA22c+dA22d+dA22e;dA(2,i)=dA22;
                    dS33 = dS33a+dS33b+dS33c+dS33d;dS(3,i)=dS33;
                    dA33 = dA33a+dA33b+dA33c+dA33d;dA(3,i)=dA33;
                    dS44 = dS44a+dS44b+dS44c;dS(4,i)=dS44;
                    dA44 = dA44a+dA44b+dA44c;dA(4,i)=dA44;
                    dS55 = dS55a+dS55b;dS(5,i)=dS55;
                    dA55 = dA55a+dA55b;dA(5,i)=dA55;
                    dS66 = dS66a;dS(6,i)=dS66;
                    dA66 = dA66a;dA(6,i)=dA66;
                end
            elseif shot(i) == shot(3)
                if rec(j) == shot(i)
                    continue
                end
                if rec(j) == shot(1)
                    g(2,1) = 0;g(3,1) = g(2,1);
                    p0 = sind(k0)*(Z_mod(2)+g(2,1)*A_mod(2));
                    y0 = [p0*(1/Z_mod(1))*(eps-A_mod(1)) 0 p0*(1/Z_mod(1))*A_mod(1) -1 p0*(1/Z_mod(1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(1,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(3)-rec(2))/(rec(3)-rec(1)))*RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2))+...
                        ((rec(2)-rec(1))/(rec(3)-rec(1)))*RP(j,i)*(Z_mod(1)+g(1,1)*A_mod(1));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(3)-rec(2))^2+(rec(2)-rec(1))^2)*(RP(j,i)/(rec(3)-rec(1)))^2;
                    dS222a = dt*((rec(3)-rec(2))/(rec(3)-rec(1)))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA222a = dt*((rec(3)-rec(2))/(rec(3)-rec(1)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS111a = dt*((rec(2)-rec(1))/(rec(3)-rec(1)))*RP(j,i)/((1+g(1,1)^2)*dd);
                    dA111a = dt*((rec(2)-rec(1))/(rec(3)-rec(1)))*RP(j,i)*g(1,1)/((1+g(1,1)^2)*dd); 
                elseif rec(j) == shot(2)
                    g(2,1)  = 0;
                    Tt(j,i) = RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dS222b = dt*((rec(3)-rec(2))/(rec(3)-rec(2)))*RP(j,i)/((1+g(2,1)^2)*RP(j,i)^2);
                    dA222b = dt*((rec(3)-rec(2))/(rec(3)-rec(2)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*RP(j,i)^2);            
                elseif rec(j) == shot(4)
                    g(3,1) = g(2,1);
                    Tt(j,i) = RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dS333a = dt*((rec(4)-rec(3))/(rec(4)-rec(3)))*RP(j,i)/((1+g(3,1)^2)*RP(j,i)^2);
                    dA333a = dt*((rec(4)-rec(3))/(rec(4)-rec(3)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*RP(j,i)^2);
                elseif rec(j) == shot(5)
                    p0 = sind(teta(j,i))*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(4)-rec(3))/(rec(5)-rec(3)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(5)-rec(3)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(4)-rec(3))^2+(rec(5)-rec(4))^2)*(RP(j,i)/(rec(5)-rec(3)))^2;
                    dS333b = dt*((rec(4)-rec(3))/(rec(5)-rec(3)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA333b = dt*((rec(4)-rec(3))/(rec(5)-rec(3)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS444a = dt*((rec(5)-rec(4))/(rec(5)-rec(3)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA444a = dt*((rec(5)-rec(4))/(rec(5)-rec(3)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                elseif rec(j) == shot(6)
                    p0 = sind(teta(j,i))*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(5))*(eps-A_mod(5)) 0 p0*(1/Z_mod(4))*A_mod(5) -1 p0*(1/Z_mod(5))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(5,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(4)-rec(3))/(rec(6)-rec(3)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(6)-rec(3)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(6)-rec(3)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(4)-rec(3))^2+(rec(5)-rec(4))^2+(rec(6)-rec(5))^2)*(RP(j,i)/(rec(6)-rec(3)))^2;
                    dS333c = dt*((rec(4)-rec(3))/(rec(6)-rec(3)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA333c = dt*((rec(4)-rec(3))/(rec(6)-rec(3)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS444b = dt*((rec(5)-rec(4))/(rec(6)-rec(3)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA444b = dt*((rec(5)-rec(4))/(rec(6)-rec(3)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS555a = dt*((rec(6)-rec(5))/(rec(6)-rec(3)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA555a = dt*((rec(6)-rec(5))/(rec(6)-rec(3)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                elseif rec(j) == shot(7)
                    p0 = sind(teta(j,i))*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(5))*(eps-A_mod(5)) 0 p0*(1/Z_mod(5))*A_mod(5) -1 p0*(1/Z_mod(5))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(5,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(5)+g(5,1)*A_mod(5));
                    y0 = [p0*(1/Z_mod(6))*(eps-A_mod(6)) 0 p0*(1/Z_mod(6))*A_mod(6) -1 p0*(1/Z_mod(6))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    Tt(j,i) = ((rec(4)-rec(3))/(rec(7)-rec(3)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(7)-rec(3)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(7)-rec(3)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5))+...
                        ((rec(7)-rec(6))/(rec(7)-rec(3)))*RP(j,i)*(Z_mod(6)+g(6,1)*A_mod(6));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(4)-rec(3))^2+(rec(5)-rec(4))^2+(rec(6)-rec(5))^2+(rec(7)-rec(6))^2)*(RP(j,i)/(rec(7)-rec(3)))^2;
                    dS333d = dt*((rec(4)-rec(3))/(rec(7)-rec(3)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA333d = dt*((rec(4)-rec(3))/(rec(7)-rec(3)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS444c = dt*((rec(5)-rec(4))/(rec(7)-rec(3)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA444c = dt*((rec(5)-rec(4))/(rec(7)-rec(3)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS555b = dt*((rec(6)-rec(5))/(rec(7)-rec(3)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA555b = dt*((rec(6)-rec(5))/(rec(7)-rec(3)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                    dS666a = dt*((rec(7)-rec(6))/(rec(7)-rec(3)))*RP(j,i)/((1+g(6,1)^2)*dd);
                    dA666a = dt*((rec(7)-rec(6))/(rec(7)-rec(3)))*RP(j,i)*g(6,1)/((1+g(6,1)^2)*dd);
                    %RATA SIRT TEMBAKAN
                    dS111 = dS111a;dS(3,i) = dS111;
                    dA111 = dA111a;dA(3,i) = dA111;
                    dS222 = dS222a+dS222b;dS(2,i) = dS222;
                    dA222 = dA222a+dA222b;dS(2,i) = dA222;
                    dS333 = dS333a+dS333b+dS333c+dS333d;dS(3,i) = dS333;
                    dA333 = dA333a+dA333b+dA333c+dA333d;dA(3,i) = dS333;
                    dS444 = dS444a+dS444b+dS444c;dS(4,i) = dS444;
                    dA444 = dA444a+dA444b+dA444c;dA(4,i) = dA444;
                    dS555 = dS555a+dS555b;dS(5,i) = dS555;
                    dA555 = dA555a+dA555b;dA(5,i) = dA555;
                    dS666 = dS555a;dS(6,i) = dS666;
                    dA666 = dA666a;dA(6,i) = dA666;
                end
            elseif shot(i) == shot(4)
                if rec(j) == shot(i)
                    continue
                end
                if rec(j) == shot(1)
                    g(3,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(2))*(eps-A_mod(2)) 0 p0*(1/Z_mod(2))*A_mod(2) -1 p0*(1/Z_mod(2))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(2,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(2)+g(2,1)*A_mod(2));
                    y0 = [p0*(1/Z_mod(1))*(eps-A_mod(1)) 0 p0*(1/Z_mod(1))*A_mod(1) -1 p0*(1/Z_mod(1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(1,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(2)-rec(1))/(rec(4)-rec(1)))*RP(j,i)*(Z_mod(1)+g(1,1)*A_mod(1))+...
                        ((rec(3)-rec(2))/(rec(4)-rec(1)))*RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2))+...
                        ((rec(4)-rec(3))/(rec(4)-rec(1)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(2)-rec(1))^2+(rec(3)-rec(2))^2+(rec(4)-rec(3))^2)*(RP(j,i)/(rec(4)-rec(1)))^2;
                    dS411a = dt*((rec(2)-rec(1))/(rec(4)-rec(1)))*RP(j,i)/((1+g(1,1)^2)*dd);
                    dA411a = dt*((rec(2)-rec(1))/(rec(4)-rec(1)))*RP(j,i)*g(1,1)/((1+g(1,1)^2)*dd);
                    dS422a = dt*((rec(3)-rec(2))/(rec(4)-rec(1)))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA422a = dt*((rec(3)-rec(2))/(rec(4)-rec(1)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS433a = dt*((rec(4)-rec(3))/(rec(4)-rec(1)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA433a = dt*((rec(4)-rec(3))/(rec(4)-rec(1)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                elseif rec(j) == shot(2)
                    g(3,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(2))*(eps-A_mod(2)) 0 p0*(1/Z_mod(2))*A_mod(2) -1 p0*(1/Z_mod(2))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(2,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(3)-rec(2))/(rec(4)-rec(2)))*RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2))+...
                        ((rec(4)-rec(3))/(rec(4)-rec(2)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(3)-rec(2))^2+(rec(4)-rec(3))^2)*(RP(j,i)/(rec(4)-rec(2)))^2;
                    dS422b = dt*((rec(3)-rec(2))/(rec(4)-rec(2)))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA422b = dt*((rec(3)-rec(2))/(rec(4)-rec(2)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS433b = dt*((rec(4)-rec(3))/(rec(4)-rec(2)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA433b = dt*((rec(4)-rec(3))/(rec(4)-rec(2)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                elseif rec(j) == shot(3)
                    Tt(j,i) = ((rec(4)-rec(3))/(rec(4)-rec(2)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(4)-rec(3))^2)*(RP(j,i)/(rec(4)-rec(2)))^2;
                    dS433c = dt*((rec(4)-rec(3))/(rec(4)-rec(3)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA433c = dt*((rec(4)-rec(3))/(rec(4)-rec(3)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);                 
                elseif rec(j) == shot(5)
                    g(4,1) = g(3,1);
                    Tt(j,i) = ((rec(5)-rec(4))/(rec(5)-rec(4)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(5)-rec(4))^2)*(RP(j,i)/(rec(5)-rec(4)))^2;
                    DS444a = dt*((rec(5)-rec(4))/(rec(5)-rec(4)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    DA444a = dt*((rec(5)-rec(4))/(rec(5)-rec(4)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);              
                elseif rec(j) == shot(6)
                    g(4,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(5))*(eps-A_mod(5)) 0 p0*(1/Z_mod(5))*A_mod(5) -1 p0*(1/Z_mod(5))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(5,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(5)-rec(4))/(rec(6)-rec(4)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(6)-rec(4)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(5)-rec(4))^2+(rec(6)-rec(5))^2)*(RP(j,i)/(rec(6)-rec(4)))^2;
                    DS444b = dt*((rec(5)-rec(4))/(rec(6)-rec(4)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    DA444b = dt*((rec(5)-rec(4))/(rec(6)-rec(4)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS455a = dt*((rec(6)-rec(5))/(rec(6)-rec(4)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA455a = dt*((rec(6)-rec(5))/(rec(6)-rec(4)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                elseif rec(j) == shot(7)
                    g(4,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(5))*(eps-A_mod(5)) 0 p0*(1/Z_mod(5))*A_mod(5) -1 p0*(1/Z_mod(5))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(5,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(teta(j,i))*(Z_mod(5)+g(5,1)*A_mod(5));
                    y0 = [p0*(1/Z_mod(5))*(eps-A_mod(5)) 0 p0*(1/Z_mod(5))*A_mod(5) -1 p0*(1/Z_mod(5))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(6,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(5)-rec(4))/(rec(7)-rec(4)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(7)-rec(4)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5))+...
                        ((rec(7)-rec(6))/(rec(7)-rec(4)))*RP(j,i)*(Z_mod(6)+g(6,1)*A_mod(6));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(5)-rec(4))^2+(rec(6)-rec(5))^2+(rec(7)-rec(6))^2)*(RP(j,i)/(rec(7)-rec(4)))^2;
                    DS444c = dt*((rec(5)-rec(4))/(rec(7)-rec(4)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    DA444c = dt*((rec(5)-rec(4))/(rec(7)-rec(4)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS455b = dt*((rec(6)-rec(5))/(rec(7)-rec(4)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA455b = dt*((rec(6)-rec(5))/(rec(7)-rec(4)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                    dS466a = dt*((rec(7)-rec(6))/(rec(7)-rec(4)))*RP(j,i)/((1+g(6,1)^2)*dd);
                    dA466a = dt*((rec(7)-rec(6))/(rec(7)-rec(4)))*RP(j,i)*g(6,1)/((1+g(6,1)^2)*dd);
                    dS4111 = dS411a;dS(4,i)=dS4111;
                    dA4111 = dA411a;dA(4,i)=dA4111;
                    dS4222 = dS422a+dS422b;dS(2,i)=dS4222;
                    dA4222 = dA422a+dA422b;dA(2,i)=dS4222;
                    dS4333 = dS433a+dS433b+dS433c;dS(3,i)=dS4333;
                    dA4333 = dA433a+dA433b+dA433c;dA(3,i)=dS4333;
                    DS4444 = DS444a+dS444b+dS444c;dS(4,i)=DS4444;
                    DA4444 = DA444a+dA444b+dA444c;dA(4,i)=DS4444;
                    dS4555 = dS455a+dS455b;dS(5,i)=dS4555;
                    dA4555 = dA455a+dA455b;dA(5,i)=dS4555;
                    dS4666 = dS466a;dS(6,i)=dS4666;
                    dA4666 = dA466a;dA(6,i)=dA4666;
                end
            elseif shot(i) == shot(5)
                if rec(j) == shot(i)
                    continue
                end
                if rec(j) == shot(1)
                    g(4,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(3))*(eps-A_mod(3)) 0 p0*(1/Z_mod(3))*A_mod(3) -1 p0*(1/Z_mod(3))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(2))*(eps-A_mod(2)) 0 p0*(1/Z_mod(2))*A_mod(2) -1 p0*(1/Z_mod(2))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(2,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(2)+g(2,1)*A_mod(2));
                    y0 = [p0*(1/Z_mod(1))*(eps-A_mod(1)) 0 p0*(1/Z_mod(1))*A_mod(1) -1 p0*(1/Z_mod(1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(1,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(2)-rec(1))/(rec(5)-rec(1)))*RP(j,i)*(Z_mod(1)+g(1,1)*A_mod(1))+...
                        ((rec(3)-rec(2))/(rec(5)-rec(1)))*RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2))+...
                        ((rec(4)-rec(3))/(rec(5)-rec(1)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(5)-rec(1)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(5)-rec(4))^2+(rec(4)-rec(3))^2+(rec(3)-rec(2))^2+(rec(2)-rec(1))^2)*(RP(j,i)/(rec(5)-rec(1)))^2;
                    dS511a = dt*((rec(2)-rec(1))/(rec(5)-rec(1)))*RP(j,i)/((1+g(1,1)^2)*dd);
                    dA511a = dt*((rec(2)-rec(1))/(rec(5)-rec(1)))*RP(j,i)*g(1,1)/((1+g(1,1)^2)*dd);
                    dS522a = dt*((rec(3)-rec(2))/(rec(5)-rec(1)))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA522a = dt*((rec(3)-rec(2))/(rec(5)-rec(1)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS533a = dt*((rec(4)-rec(3))/(rec(5)-rec(1)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA533a = dt*((rec(4)-rec(3))/(rec(5)-rec(1)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS544a = dt*((rec(5)-rec(4))/(rec(5)-rec(1)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA544a = dt*((rec(5)-rec(4))/(rec(5)-rec(1)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                elseif rec(j) == shot(2)
                    g(4,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(3))*(eps-A_mod(3)) 0 p0*(1/Z_mod(3))*A_mod(3) -1 p0*(1/Z_mod(3))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(2))*(eps-A_mod(2)) 0 p0*(1/Z_mod(2))*A_mod(2) -1 p0*(1/Z_mod(2))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(2,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(3)-rec(2))/(rec(5)-rec(2)))*RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2))+...
                        ((rec(4)-rec(3))/(rec(5)-rec(2)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(5)-rec(2)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(5)-rec(4))^2+(rec(4)-rec(3))^2+(rec(3)-rec(2))^2)*(RP(j,i)/(rec(5)-rec(2)))^2;
                    dS522b = dt*((rec(3)-rec(2))/(rec(5)-rec(2)))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA522b = dt*((rec(3)-rec(2))/(rec(5)-rec(2)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS533b = dt*((rec(4)-rec(3))/(rec(5)-rec(2)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA533b = dt*((rec(4)-rec(3))/(rec(5)-rec(2)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS544b = dt*((rec(5)-rec(4))/(rec(5)-rec(2)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA544b = dt*((rec(5)-rec(4))/(rec(5)-rec(2)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                elseif rec(j) == shot(3)
                    g(4,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(3))*(eps-A_mod(3)) 0 p0*(1/Z_mod(3))*A_mod(3) -1 p0*(1/Z_mod(3))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(4)-rec(3))/(rec(5)-rec(3)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(5)-rec(3)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(5)-rec(4))^2+(rec(4)-rec(3))^2)*(RP(j,i)/(rec(5)-rec(3)))^2;
                    dS533c = dt*((rec(4)-rec(3))/(rec(5)-rec(3)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA533c = dt*((rec(4)-rec(3))/(rec(5)-rec(3)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS544c = dt*((rec(5)-rec(4))/(rec(5)-rec(3)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA544c = dt*((rec(5)-rec(4))/(rec(5)-rec(3)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                elseif rec(j) == shot(4)
                    g(4,1) = 0;
                    Tt(j,i) = ((rec(5)-rec(4))/(rec(5)-rec(4)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(5)-rec(4))^2)*(RP(j,i)/(rec(5)-rec(4)))^2;
                    dS544d = dt*((rec(5)-rec(4))/(rec(5)-rec(4)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA544d = dt*((rec(5)-rec(4))/(rec(5)-rec(4)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                elseif rec(j) == shot(6)
                    g(5,1) = g(4,1);
                    Tt(j,i) = ((rec(6)-rec(5))/(rec(6)-rec(5)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(6)-rec(5))^2)*(RP(j,i)/(rec(6)-rec(5)))^2;
                    DS555a = dt*((rec(6)-rec(5))/(rec(6)-rec(5)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    DA555a = dt*((rec(6)-rec(5))/(rec(6)-rec(5)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                elseif rec(j) == shot(7)
                    g(5,1) = g(4,1);
                    p0 = sind(teta(j,i))*(Z_mod(5)+g(5,1)*A_mod(5));
                    y0 = [p0*(1/Z_mod(6))*(eps-A_mod(6)) 0 p0*(1/Z_mod(6))*A_mod(6) -1 p0*(1/Z_mod(6))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(6,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(6)-rec(5))/(rec(7)-rec(5)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5))+...
                        ((rec(7)-rec(6))/(rec(7)-rec(5)))*RP(j,i)*(Z_mod(6)+g(6,1)*A_mod(6));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(6)-rec(5))^2+(rec(7)-rec(5))^2)*(RP(j,i)/(rec(7)-rec(5)))^2;
                    DS555b = dt*((rec(6)-rec(5))/(rec(7)-rec(5)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    DA555b = dt*((rec(6)-rec(5))/(rec(7)-rec(5)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                    dS566a = dt*((rec(7)-rec(6))/(rec(7)-rec(5)))*RP(j,i)/((1+g(6,1)^2)*dd);
                    dA566a = dt*((rec(7)-rec(6))/(rec(7)-rec(5)))*RP(j,i)*g(6,1)/((1+g(6,1)^2)*dd);
                    dS511 = dS511a;dS(1,i) = dS511;
                    dA511 = dA511a;dA(1,i) = dA511;
                    dS522 = dS522a+dS522b;dS(2,i) = dS522;
                    dA522 = dA522a+dA522b;dA(2,i) = dA522;
                    dS533 = dS533a+dS533b+dS533c;dS(3,i) = dS533;
                    dA533 = dA533a+dA533b+dA533c;dA(3,i) = dA533;
                    dS544 = dS544a+dS544b+dS544c+dS544d;dS(4,i) = dS544;
                    dA544 = dA544a+dA544b+dA544c+dA544d;dA(4,i) = dA544;
                    DS555 = dS55a+dS55b;dS(5,i) = DS555;
                    DA555 = dA55a+dA55b;dA(5,i) = DA555;
                    dS566 = dS566a;dS(6,i) = dS566;
                    dA566 = dA566a;dA(6,i) = dA566;
                end
            elseif shot(i) == shot(6)
                if rec(j) == shot(i)
                    continue
                end
                if rec(j) == shot(1)
                    g(5,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(5)+g(5,1)*A_mod(5));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(3))*(eps-A_mod(3)) 0 p0*(1/Z_mod(3))*A_mod(3) -1 p0*(1/Z_mod(3))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(2))*(eps-A_mod(2)) 0 p0*(1/Z_mod(2))*A_mod(2) -1 p0*(1/Z_mod(2))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(2,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(2)+g(2,1)*A_mod(2));
                    y0 = [p0*(1/Z_mod(1))*(eps-A_mod(1)) 0 p0*(1/Z_mod(1))*A_mod(1) -1 p0*(1/Z_mod(1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(1,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(2)-rec(1))/(rec(6)-rec(1)))*RP(j,i)*(Z_mod(1)+g(1,1)*A_mod(1))+...
                        ((rec(3)-rec(2))/(rec(6)-rec(1)))*RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2))+...
                        ((rec(4)-rec(3))/(rec(6)-rec(1)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(6)-rec(1)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(6)-rec(1)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(2)-rec(1))^2+(rec(3)-rec(2))^2+(rec(4)-rec(3))^2+(rec(5)-rec(4))^2+(rec(6)-rec(5))^2)*(RP(j,i)/(rec(6)-rec(1)))^2;
                    dS611a = dt*((rec(2)-rec(1))/(rec(6)-rec(1)))*RP(j,i)/((1+g(1,1)^2)*dd);
                    dA611a = dt*((rec(2)-rec(1))/(rec(6)-rec(1)))*RP(j,i)*g(1,1)/((1+g(1,1)^2)*dd);
                    dS622a = dt*((rec(3)-rec(2))/(rec(6)-rec(1)))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA622a = dt*((rec(3)-rec(2))/(rec(6)-rec(1)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS633a = dt*((rec(4)-rec(3))/(rec(6)-rec(1)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA633a = dt*((rec(4)-rec(3))/(rec(6)-rec(1)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS644a = dt*((rec(5)-rec(4))/(rec(6)-rec(1)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA644a = dt*((rec(5)-rec(4))/(rec(6)-rec(1)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS655a = dt*((rec(6)-rec(5))/(rec(6)-rec(1)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA655a = dt*((rec(6)-rec(5))/(rec(6)-rec(1)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                elseif rec(j) == shot(2)
                    g(5,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(5)+g(5,1)*A_mod(5));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(3))*(eps-A_mod(3)) 0 p0*(1/Z_mod(3))*A_mod(3) -1 p0*(1/Z_mod(3))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(2))*(eps-A_mod(2)) 0 p0*(1/Z_mod(2))*A_mod(2) -1 p0*(1/Z_mod(2))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(2,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(3)-rec(2))/(rec(6)-rec(2)))*RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2))+...
                        ((rec(4)-rec(3))/(rec(6)-rec(2)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(6)-rec(2)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(6)-rec(2)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(3)-rec(2))^2+(rec(4)-rec(3))^2+(rec(5)-rec(4))^2+(rec(6)-rec(5))^2)*(RP(j,i)/(rec(6)-rec(2)))^2;
                    dS622b = dt*((rec(3)-rec(2))/(rec(6)-rec(2)))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA622b = dt*((rec(3)-rec(2))/(rec(6)-rec(2)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS633b = dt*((rec(4)-rec(3))/(rec(6)-rec(2)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA633b = dt*((rec(4)-rec(3))/(rec(6)-rec(2)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS644b = dt*((rec(5)-rec(4))/(rec(6)-rec(2)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA644b = dt*((rec(5)-rec(4))/(rec(6)-rec(2)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS655b = dt*((rec(6)-rec(5))/(rec(6)-rec(2)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA655b = dt*((rec(6)-rec(5))/(rec(6)-rec(2)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                elseif rec(j) == shot(3)
                    g(5,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(5)+g(5,1)*A_mod(5));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(3))*(eps-A_mod(3)) 0 p0*(1/Z_mod(3))*A_mod(3) -1 p0*(1/Z_mod(3))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(4)-rec(3))/(rec(6)-rec(3)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(6)-rec(3)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(6)-rec(3)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(4)-rec(3))^2+(rec(5)-rec(4))^2+(rec(6)-rec(5))^2)*(RP(j,i)/(rec(6)-rec(3)))^2;
                    dS633c = dt*((rec(4)-rec(3))/(rec(6)-rec(3)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA633c = dt*((rec(4)-rec(3))/(rec(6)-rec(3)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS644c = dt*((rec(5)-rec(4))/(rec(6)-rec(3)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA644c = dt*((rec(5)-rec(4))/(rec(6)-rec(3)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS655c = dt*((rec(6)-rec(5))/(rec(6)-rec(3)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA655c = dt*((rec(6)-rec(5))/(rec(6)-rec(3)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                elseif rec(j) == shot(4)
                    g(5,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(5)+g(5,1)*A_mod(5));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(5)-rec(4))/(rec(6)-rec(4)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(6)-rec(4)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(5)-rec(4))^2+(rec(6)-rec(5))^2)*(RP(j,i)/(rec(6)-rec(4)))^2;
                    dS644d = dt*((rec(5)-rec(4))/(rec(6)-rec(4)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA644d = dt*((rec(5)-rec(4))/(rec(6)-rec(4)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS655d = dt*((rec(6)-rec(5))/(rec(6)-rec(4)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA655d = dt*((rec(6)-rec(5))/(rec(6)-rec(4)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                elseif rec(j) == shot(5)
                    g(5,1) = 0;
                    Tt(j,i) = ((rec(6)-rec(5))/(rec(6)-rec(5)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(6)-rec(5))^2)*(RP(j,i)/(rec(6)-rec(5)))^2;
                    dS655e = dt*((rec(6)-rec(5))/(rec(6)-rec(5)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA655e = dt*((rec(6)-rec(5))/(rec(6)-rec(5)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                elseif rec(j) == shot(7)
                    g(6,1) =g(5,1);
                    Tt(j,i) = ((rec(7)-rec(6))/(rec(7)-rec(6)))*RP(j,i)*(Z_mod(6)+g(6,1)*A_mod(6));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(7)-rec(6))^2)*(RP(j,i)/(rec(7)-rec(6)))^2;
                    DS666a = dt*((rec(7)-rec(6))/(rec(7)-rec(6)))*RP(j,i)/((1+g(6,1)^2)*dd);
                    DA666a = dt*((rec(7)-rec(6))/(rec(7)-rec(6)))*RP(j,i)*g(6,1)/((1+g(6,1)^2)*dd);
                    dS611 = dS611a;dS(1,i) = dS611;
                    dA611 = dA611a;dA(1,i) = dA611;
                    dS622 = dS622a+dS622b;dS(2,i) = dS622;
                    dA622 = dA622a+dA622b;dA(2,i) = dA622;
                    dS633 = dS633a+dS633b+dS633c;dS(3,i) = dS633;
                    dA633 = dA633a+dA633b+dA633c;dA(3,i) = dA633;
                    dS644 = dS644a+dS644b+dS644c+dS644d;dS(4,i) = dS644;
                    dA644 = dA644a+dA644b+dA644c+dA644d;dA(4,i) = dA644;
                    dS655 = dS655a+dS655b+dS655c+dS655d+dS655e;dS(5,i) = dS655;
                    dA655 = dA655a+dA655b+dA655c+dA655d+dA655e;dA(5,i) = dA655;
                    DS666 = DS666a;dS(6,i) = DS666;
                    DA666 = DA666a;dA(6,i) = DA666;
                end
            elseif shot(i) == shot(7)
                if rec(j) == shot(i)
                    continue
                end
                if rec(j) == shot(1)
                    g(6,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(6)+g(6,1)*A_mod(6));
                    y0 = [p0*(1/Z_mod(5))*(eps-A_mod(5)) 0 p0*(1/Z_mod(5))*A_mod(5) -1 p0*(1/Z_mod(5))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(5,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(5)+g(5,1)*A_mod(5));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(3))*(eps-A_mod(3)) 0 p0*(1/Z_mod(3))*A_mod(3) -1 p0*(1/Z_mod(3))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(2))*(eps-A_mod(2)) 0 p0*(1/Z_mod(2))*A_mod(2) -1 p0*(1/Z_mod(2))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(2,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(2)+g(2,1)*A_mod(2));
                    y0 = [p0*(1/Z_mod(1))*(eps-A_mod(1)) 0 p0*(1/Z_mod(1))*A_mod(1) -1 p0*(1/Z_mod(1))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(1,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(2)-rec(1))/(rec(7)-rec(1)))*RP(j,i)*(Z_mod(1)+g(1,1)*A_mod(1))+...
                        ((rec(3)-rec(2))/(rec(7)-rec(1)))*RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2))+...
                        ((rec(4)-rec(3))/(rec(7)-rec(1)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(7)-rec(1)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(7)-rec(1)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5))+...
                        ((rec(7)-rec(6))/(rec(7)-rec(1)))*RP(j,i)*(Z_mod(6)+g(6,1)*A_mod(6));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(7)-rec(6))^2+(rec(6)-rec(5))^2+(rec(5)-rec(4))^2+(rec(4)-rec(3))^2+(rec(3)-rec(2))^2+(rec(2)-rec(1))^2)*(RP(j,i)/(rec(7)-rec(1)))^2;
                    dS711a = dt*((rec(2)-rec(1))/(rec(7)-rec(1)))*RP(j,i)/((1+g(1,1)^2)*dd);
                    dA711a = dt*((rec(2)-rec(1))/(rec(7)-rec(1)))*RP(j,i)*g(1,1)/((1+g(1,1)^2)*dd);
                    dS722a = dt*((rec(3)-rec(2))/(rec(7)-rec(1)))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA722a = dt*((rec(3)-rec(2))/(rec(7)-rec(1)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS733a = dt*((rec(4)-rec(3))/(rec(7)-rec(1)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA733a = dt*((rec(4)-rec(3))/(rec(7)-rec(1)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS744a = dt*((rec(5)-rec(4))/(rec(7)-rec(1)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA744a = dt*((rec(5)-rec(4))/(rec(7)-rec(1)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS755a = dt*((rec(6)-rec(5))/(rec(7)-rec(1)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA755a = dt*((rec(6)-rec(5))/(rec(7)-rec(1)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                    dS766a = dt*((rec(7)-rec(6))/(rec(7)-rec(1)))*RP(j,i)/((1+g(6,1)^2)*dd);
                    dA766a = dt*((rec(7)-rec(6))/(rec(7)-rec(1)))*RP(j,i)*g(6,1)/((1+g(6,1)^2)*dd);
                elseif rec(j) == shot(2)
                    g(6,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(6)+g(6,1)*A_mod(6));
                    y0 = [p0*(1/Z_mod(5))*(eps-A_mod(5)) 0 p0*(1/Z_mod(5))*A_mod(5) -1 p0*(1/Z_mod(5))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(5,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(5)+g(5,1)*A_mod(5));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(3))*(eps-A_mod(3)) 0 p0*(1/Z_mod(3))*A_mod(3) -1 p0*(1/Z_mod(3))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(3)+g(3,1)*A_mod(3));
                    y0 = [p0*(1/Z_mod(2))*(eps-A_mod(2)) 0 p0*(1/Z_mod(2))*A_mod(2) -1 p0*(1/Z_mod(2))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(2,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(3)-rec(2))/(rec(7)-rec(2)))*RP(j,i)*(Z_mod(2)+g(2,1)*A_mod(2))+...
                        ((rec(4)-rec(3))/(rec(7)-rec(2)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(7)-rec(2)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(7)-rec(2)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5))+...
                        ((rec(7)-rec(6))/(rec(7)-rec(2)))*RP(j,i)*(Z_mod(6)+g(6,1)*A_mod(6));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(7)-rec(6))^2+(rec(6)-rec(5))^2+(rec(5)-rec(4))^2+(rec(4)-rec(3))^2+(rec(3)-rec(2))^2)*(RP(j,i)/(rec(7)-rec(2)))^2;
                    dS722b = dt*((rec(3)-rec(2))/(rec(7)-rec(2)))*RP(j,i)/((1+g(2,1)^2)*dd);
                    dA722b = dt*((rec(3)-rec(2))/(rec(7)-rec(2)))*RP(j,i)*g(2,1)/((1+g(2,1)^2)*dd);
                    dS733b = dt*((rec(4)-rec(3))/(rec(7)-rec(2)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA733b = dt*((rec(4)-rec(3))/(rec(7)-rec(2)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS744b = dt*((rec(5)-rec(4))/(rec(7)-rec(2)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA744b = dt*((rec(5)-rec(4))/(rec(7)-rec(2)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS755b = dt*((rec(6)-rec(5))/(rec(7)-rec(2)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA755b = dt*((rec(6)-rec(5))/(rec(7)-rec(2)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                    dS766b = dt*((rec(7)-rec(6))/(rec(7)-rec(2)))*RP(j,i)/((1+g(6,1)^2)*dd);
                    dA766b = dt*((rec(7)-rec(6))/(rec(7)-rec(2)))*RP(j,i)*g(6,1)/((1+g(6,1)^2)*dd);
                elseif rec(j) == shot(3)
                    g(6,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(6)+g(6,1)*A_mod(6));
                    y0 = [p0*(1/Z_mod(5))*(eps-A_mod(5)) 0 p0*(1/Z_mod(5))*A_mod(5) -1 p0*(1/Z_mod(5))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(5,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(5)+g(5,1)*A_mod(5));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(4)+g(4,1)*A_mod(4));
                    y0 = [p0*(1/Z_mod(3))*(eps-A_mod(3)) 0 p0*(1/Z_mod(3))*A_mod(3) -1 p0*(1/Z_mod(3))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(3,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) = ((rec(4)-rec(3))/(rec(7)-rec(3)))*RP(j,i)*(Z_mod(3)+g(3,1)*A_mod(3))+...
                        ((rec(5)-rec(4))/(rec(7)-rec(3)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(7)-rec(3)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5))+...
                        ((rec(7)-rec(6))/(rec(7)-rec(3)))*RP(j,i)*(Z_mod(6)+g(6,1)*A_mod(6));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(7)-rec(6))^2+(rec(6)-rec(5))^2+(rec(5)-rec(4))^2+(rec(4)-rec(3))^2)*(RP(j,i)/(rec(7)-rec(3)))^2;
                    dS733c = dt*((rec(4)-rec(3))/(rec(7)-rec(3)))*RP(j,i)/((1+g(3,1)^2)*dd);
                    dA733c = dt*((rec(4)-rec(3))/(rec(7)-rec(3)))*RP(j,i)*g(3,1)/((1+g(3,1)^2)*dd);
                    dS744c = dt*((rec(5)-rec(4))/(rec(7)-rec(3)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA744c = dt*((rec(5)-rec(4))/(rec(7)-rec(3)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS755c = dt*((rec(6)-rec(5))/(rec(7)-rec(3)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA755c = dt*((rec(6)-rec(5))/(rec(7)-rec(3)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                    dS766c = dt*((rec(7)-rec(6))/(rec(7)-rec(3)))*RP(j,i)/((1+g(6,1)^2)*dd);
                    dA766c = dt*((rec(7)-rec(6))/(rec(7)-rec(3)))*RP(j,i)*g(6,1)/((1+g(6,1)^2)*dd);
                elseif rec(j) == shot(4)
                    g(6,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(6)+g(6,1)*A_mod(6));
                    y0 = [p0*(1/Z_mod(5))*(eps-A_mod(5)) 0 p0*(1/Z_mod(5))*A_mod(5) -1 p0*(1/Z_mod(5))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(5,1) = -1*sind(k0)^2*cosd(k0)^2;
                    p0 = sind(k0)*(Z_mod(5)+g(5,1)*A_mod(5));
                    y0 = [p0*(1/Z_mod(4))*(eps-A_mod(4)) 0 p0*(1/Z_mod(4))*A_mod(4) -1 p0*(1/Z_mod(4))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(4,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) =((rec(5)-rec(4))/(rec(7)-rec(4)))*RP(j,i)*(Z_mod(4)+g(4,1)*A_mod(4))+...
                        ((rec(6)-rec(5))/(rec(7)-rec(4)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5))+...
                        ((rec(7)-rec(6))/(rec(7)-rec(4)))*RP(j,i)*(Z_mod(6)+g(6,1)*A_mod(6));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(7)-rec(6))^2+(rec(6)-rec(5))^2+(rec(5)-rec(4))^2)*(RP(j,i)/(rec(7)-rec(4)))^2;
                    dS744d = dt*((rec(5)-rec(4))/(rec(7)-rec(4)))*RP(j,i)/((1+g(4,1)^2)*dd);
                    dA744d = dt*((rec(5)-rec(4))/(rec(7)-rec(4)))*RP(j,i)*g(4,1)/((1+g(4,1)^2)*dd);
                    dS755d = dt*((rec(6)-rec(5))/(rec(7)-rec(4)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA755d = dt*((rec(6)-rec(5))/(rec(7)-rec(4)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                    dS766d = dt*((rec(7)-rec(6))/(rec(7)-rec(4)))*RP(j,i)/((1+g(6,1)^2)*dd);
                    dA766d = dt*((rec(7)-rec(6))/(rec(7)-rec(4)))*RP(j,i)*g(6,1)/((1+g(6,1)^2)*dd);
                elseif rec(j) == shot(5)
                    g(6,1) = 0;
                    p0 = sind(teta(j,i))*(Z_mod(6)+g(6,1)*A_mod(6));
                    y0 = [p0*(1/Z_mod(5))*(eps-A_mod(5)) 0 p0*(1/Z_mod(5))*A_mod(5) -1 p0*(1/Z_mod(5))];
                    k0 = abs(asind(roots(y0))); 
                    if length(k0) > 1
                        for jj = 1:length(k0)
                            k00 = k0(jj)-teta(j,i);
                        end
                        [val,index] = min(k00);k0 = k0(index);
                    end
                    g(5,1) = -1*sind(k0)^2*cosd(k0)^2;
                    Tt(j,i) =((rec(6)-rec(5))/(rec(7)-rec(5)))*RP(j,i)*(Z_mod(5)+g(5,1)*A_mod(5))+...
                        ((rec(7)-rec(6))/(rec(7)-rec(5)))*RP(j,i)*(Z_mod(6)+g(6,1)*A_mod(6));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(7)-rec(6))^2+(rec(6)-rec(5))^2)*(RP(j,i)/(rec(7)-rec(5)))^2;
                    dS755e = dt*((rec(6)-rec(5))/(rec(7)-rec(5)))*RP(j,i)/((1+g(5,1)^2)*dd);
                    dA755e = dt*((rec(6)-rec(5))/(rec(7)-rec(5)))*RP(j,i)*g(5,1)/((1+g(5,1)^2)*dd);
                    dS766e = dt*((rec(7)-rec(6))/(rec(7)-rec(5)))*RP(j,i)/((1+g(6,1)^2)*dd);
                    dA766e = dt*((rec(7)-rec(6))/(rec(7)-rec(5)))*RP(j,i)*g(6,1)/((1+g(6,1)^2)*dd);
                elseif rec(j) == shot(6)
                    g(6,1) = 0;
                    Tt(j,i) =((rec(7)-rec(6))/(rec(7)-rec(6)))*RP(j,i)*(Z_mod(6)+g(6,1)*A_mod(6));
                    dt = TT(j,i)-Tt(j,i); D(j,i) = dt;
                    dd =((rec(7)-rec(6))^2)*(RP(j,i)/(rec(7)-rec(6)))^2;
                    dS766f = dt*((rec(7)-rec(6))/(rec(7)-rec(6)))*RP(j,i)/((1+g(6,1)^2)*dd);
                    dA766f = dt*((rec(7)-rec(6))/(rec(7)-rec(6)))*RP(j,i)*g(6,1)/((1+g(6,1)^2)*dd);
                    dS711 = dS711a;dS(1,i) = dS711;
                    dA711 = dA711a;dA(1,i) = dA711;
                    dS722 = dS722a+dS722b;dS(2,i) = dS722;
                    dA722 = dA722a+dA722b;dA(2,i) = dA722;
                    dS733 = dS733a+dS733b+dS733c;dS(3,i) = dS733;
                    dA733 = dA733a+dA733b+dA733c;dA(3,i) = dA733;
                    dS744 = dS744a+dS744b+dS744c+dS744d;dS(4,i) = dS744;
                    dA744 = dA744a+dA744b+dA744c+dA744d;dA(4,i) = dA744;
                    dS755 = dS755a+dS755b+dS755c+dS755d+dS755e;dS(5,i) = dS755;
                    dA755 = dA755a+dA755b+dA755c+dA755d+dA755e;dA(5,i) = dA755;
                    dS766 = dS766a+dS766b+dS766c+dS766d+dS766e+dS766f;dS(6,i) = dS766;
                    dA766 = dA766a+dA766b+dA766c+dA766d+dA766e+dA766f;dA(6,i) = dA766;
                end
            end
        end
         PERTU_S(1) = sum(dS(1,:))/12;PERTU_S(2)=sum(dS(2,:))/15;PERTU_S(3)=sum(dS(3,:))/25;PERTU_S(4)=sum(dS(4,:))/24;PERTU_S(5)=sum(dS(5,:))/20;PERTU_S(6)=sum(dS(6,:))/12;
         PERTU_A(1) = sum(dA(1,:))/12;PERTU_A(2)=sum(dA(2,:))/15;PERTU_A(3)=sum(dA(3,:))/25;PERTU_A(4)=sum(dA(4,:))/24;PERTU_A(5)=sum(dA(5,:))/20;PERTU_A(6)=sum(dA(6,:))/12;
        Z_mod = Z_mod+PERTU_S;A_mod = A_mod+PERTU_A;
    end
    p = p+1;
    dt = TT-Tt;
    ERROR=(mean(sum((dt).^2)))^0.5;
    if ERROR < Toleransi
        break
    end
    if ERROR > TOLERANSI_UP
        break
    end
end
n_hor = 30;
fprintf('ITERASI INVERSI SIRT : %d\n',p)

figure(2)
S = Z_mod;S(7) = [];V = 1./S;V_mod=repmat(V,1,n_hor);
V_add =V_mod(6,:);V_add_mod =vertcat(V_mod,V_add); 
Y_grid = shot;
[X,Y]= meshgrid(linspace(0,Offset,n_hor),Y_grid);
pcolor(X,Y,V_add_mod);hold on
colorbar;
set(gca,'Ydir','reverse');
title('Model Inversi Kecepatan SIRT');
xlabel('Offset Shot(m)');
ylabel('Depth Receiver(m)');
for p = 1:length(shot)
    plot(0,shot(p),'v','MarkerSize',5,'MarkerFaceColor','r');hold on
end
for p = 1:length(rec)
    plot(Offset,rec(p),'O','MarkerSize',8,'MarkerFaceColor','r');hold on
end
hold on

figure(3)
S = A_mod;S(7) = [];V = 1./S;V_mod=repmat(V,1,n_hor);
V_add =V_mod(6,:);V_add_mod =vertcat(V_mod,V_add); 
Y_grid = shot;
[X,Y]= meshgrid(linspace(0,Offset,n_hor),Y_grid);
pcolor(X,Y,V_add_mod);hold on
colorbar;
set(gca,'Ydir','reverse');
title('Delta Anisotropi SIRT');
xlabel('Offset Shot(m)');
ylabel('Depth Receiver(m)');
for p = 1:length(shot)
    plot(0,shot(p),'v','MarkerSize',5,'MarkerFaceColor','r');hold on
end
for p = 1:length(rec)
    plot(Offset,rec(p),'O','MarkerSize',8,'MarkerFaceColor','r');hold on
end
hold on




