%GEOTOMOGRAFI ART INVERSION
%KELOMPOK 1 UAS 2019/2020 
%FIXED MODEL (4 layers and 6 receivers DOWNHOLE ANISOTROPI)
%SHORT INFO---->VARIASI HANYA PADA JUMLAH SHOT,OFFSET,KECEPATAN LAPISAN,DAN
%KETEBALAN (Dengan syarat tebal lapisan pertama tidak boleh kurang dari
%posisi receiver 1 dan jumlah semua tebal lapisan sama dengan Depth)
clear;clc;close all
Depth = 30;                  
Offset = 30;            
R = linspace(5,Depth,6);            % Receiver Offset
n_shot = 20;                        % JUMLAH SHOT
shot = linspace(0,Offset,n_shot);   % Shot Offset
tbl = [10,10,5,5];                  % Tebal per lapisan
w  = 0;                             % Random Noise Ratio
dlap = [0,0,1,2];                   % Delta Anisotropi
eps = 2;                            % Epsilon Anisotropi
v = [2500,500,1000,400];slap = 1./v;
bounds = [tbl(1) tbl(1)+tbl(2) tbl(1)+tbl(2)+tbl(3) sum(tbl)];     %Batas Lapisan
if dlap(1) ~= 0
    fprintf('LAPISAN 1 ANISOTROPI\n')
else
    fprintf('LAPISAN 1 ISOTROPI\n')
end
if dlap(2) ~= 0
    fprintf('LAPISAN 2 ANISOTROPI\n')
else
    fprintf('LAPISAN 2 ISOTROPI\n')
end
if dlap(3) ~= 0
    fprintf('LAPISAN 3 ANISOTROPI\n')
else
    fprintf('LAPISAN 3 ISOTROPI\n')
end
if dlap(4) ~= 0
    fprintf('LAPISAN 4 ANISOTROPI\n')
else
    fprintf('LAPISAN 4 ISOTROPI\n')
end
if length(tbl)~= 4
    error('Jumlah lapisan harus empat\n')
end
if max(bounds) ~= Depth
    error('Jumlah tebal lapisan tidak sama dengan depth model fix\n')
end
if min(bounds)<min(R)
    error('Receiver 1 harus berada sebelum/tepat di batas lapisan 1(bounds 1)\n')
end
fprintf('VARIASI PADA SIFAT ANISO LAPISAN 1 TIDAK BERPENGARUH ATAU SAMA SAJA KARENA SUDUT TEMBAKAN SHOT PERMUKAAN\nDIASUMSIKAN SEBAGAI SUDUT FASE ANISOTROPI atau ISOTROPI PADA LAPISAN PERTAMA.\n')
%RAYPATH
RR = zeros(length(R),length(shot));teta = zeros(size(RR));
for i = 1:length(shot)
    for j = 1:length(R)
        RR(j,i) = ((shot(i))^2+(R(j))^2)^0.5;
        teta(j,i) = 90-acosd(shot(i)-RR(j,i));
    end
end
%TRAVEL TIMES
TT = zeros(size(RR));
fprintf('TRAVEL TIMES MODEL FORWARD (PURE)-->TTA\nRAYPATH LENGTH MODEL FORWARD-->RR\n')
fprintf('RANDOM NOISE RATIO : %d persen\n',w)
fprintf('TRAVEL TIMES MODEL FORWARD (NOISE)-->TT\n')
for i = 1:length(shot)
    for j = 1:length(R)
        if R(j)<=bounds(1)
            TT(j,i) = RR(j,i)/v(1);
        elseif R(j)<=bounds(2)
            if dlap(1) == 0 && dlap(2) == 0 || dlap(1)~=0 && dlap(2)==0
                TT(j,i) = ((R(j)-bounds(1))/R(j))*RR(j,i)/v(2) + (bounds(1)/R(j))*RR(j,i)/v(1);
            elseif dlap(1)==0 && dlap(2)~=0 || dlap(1)~=0 && dlap(2)~=0 
                p2 = sind(teta(j,i))*slap(1);
                y2 = [p2*v(2)*(eps-dlap(2)) 0 p2*v(2)*dlap(2) -1 p2*v(2)];
                k2 = abs(asind(roots(y2)));
                if length(k2) > 1
                    for jj = 1:length(k2)                                       % DARI 1 (Iso/Aniso) ke 2 (Aniso)
                        k22 = k2(jj)-teta(j,i);
                    end
                        [val,index] = min(k22);k2 = k2(index);
                end
                TT(j,i) = ((R(j)-bounds(1))/R(j))*RR(j,i)*(1+(-1*sind(k2)^2*cosd(k2)^2*dlap(2)))*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1);
            end
        elseif R(j)<=bounds(3)
            if dlap(1)==0 && dlap(2)==0 && dlap(3)==0 || dlap(1)~=0 && dlap(2)==0 && dlap(3)==0
                TT(j,i) = ((bounds(2)-bounds(1))/R(j))*RR(j,i)*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1)+...
                    ((R(j)-bounds(2))/R(j))*RR(j,i)*slap(3);
            elseif dlap(1)==0 && dlap(2)~=0 && dlap(3)==0 || dlap(1)~=0 && dlap(2)~=0 && dlap(3)==0
                TT(j,i) = ((bounds(2)-bounds(1))/R(j))*RR(j,i)*(1+(-1*sind(k2)^2*cosd(k2)^2*dlap(2)))*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1)+...
                    ((R(j)-bounds(2))/R(j))*RR(j,i)*slap(3);
            elseif dlap(1)==0 && dlap(2)==0 && dlap(3)~=0 || dlap(1)~=0 && dlap(2)==0 && dlap(3)~=0
                p3 = sind(teta(j,i))*slap(2);
                y3 = [p3*v(3)*(eps-dlap(3)) 0 p3*v(3)*dlap(3) -1 p3*v(3)];          % DARI 2 (Iso) ke 3 (Aniso)
                k3 = abs(asind(roots(y3)));
                if length(k3) > 1
                    for jj = 1:length(k3)
                        k33 = k3(jj)-teta(j,i);
                    end
                        [val,index] = min(k33);k3 = k3(index);
                end
                TT(j,i) = ((bounds(2)-bounds(1))/R(j))*RR(j,i)*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1)+...
                    ((R(j)-bounds(2))/R(j))*RR(j,i)*(1+(-1*sind(k3)^2*cosd(k3)^2*dlap(3)))*slap(3);
            elseif dlap(1)==0 && dlap(2)~=0 && dlap(3)~=0 || dlap(1)~=0 && dlap(2)~=0 && dlap(3)~=0
                p23 = sind(k2)*slap(2)*(1+(-1*sind(k2)^2*cosd(k2)^2*dlap(2)));
                y23 = [p23*v(3)*(eps-dlap(3)) 0 p23*v(3)*dlap(3) -1 p23*v(3)];         % DARI 2 (Aniso) ke 3 (Aniso)
                k23 = abs(asind(roots(y23)));
                if length(k23) > 1
                    for jj = 1:length(k23)
                        k223 = k23(jj)-k2;
                    end
                        [val,index] = min(k223);k23 = k23(index);
                end
                 TT(j,i) = ((bounds(2)-bounds(1))/R(j))*RR(j,i)*(1+(-1*sind(k2)^2*cosd(k2)^2*dlap(2)))*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1)+...
                   ((R(j)-bounds(2))/R(j))*RR(j,i)*(1+(-1*sind(k23)^2*cosd(k23)^2*dlap(3)))*slap(3);
            end
        elseif R(j)<=bounds(4)
            if dlap(1)==0 && dlap(2)==0 && dlap(3)==0 && dlap(4)==0 || dlap(1)~=0 && dlap(2)==0 && dlap(3)==0 && dlap(4)==0
                 TT(j,i) = ((bounds(2)-bounds(1))/R(j))*RR(j,i)*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1)+...
                    ((bounds(3)-bounds(2))/R(j))*RR(j,i)*slap(3) + ((R(j)-bounds(3))/R(j))*RR(j,i)*slap(4);
            elseif dlap(1)==0 && dlap(2)~=0 && dlap(3)==0 && dlap(4)==0 || dlap(1)~=0 && dlap(2)~=0 && dlap(3)==0 && dlap(4)==0
                 TT(j,i) = ((bounds(2)-bounds(1))/R(j))*RR(j,i)*(1+(-1*sind(k2)^2*cosd(k2)^2*dlap(2)))*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1)+...
                    ((bounds(3)-bounds(2))/R(j))*RR(j,i)*slap(3) + ((R(j)-bounds(3))/R(j))*RR(j,i)*slap(4);
            elseif dlap(1)==0 && dlap(2)==0 && dlap(3)~=0 && dlap(4)==0 || dlap(1)~=0 && dlap(2)==0 && dlap(3)~=0 && dlap(4)==0
                TT(j,i) = ((bounds(2)-bounds(1))/R(j))*RR(j,i)*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1)+...
                    ((bounds(3)-bounds(2))/R(j))*RR(j,i)*(1+(-1*sind(k3)^2*cosd(k3)^2*dlap(3)))*slap(3) + ((R(j)-bounds(3))/R(j))*RR(j,i)*slap(4);
            elseif dlap(1)==0 && dlap(2)==0 && dlap(3)==0 && dlap(4)~=0 || dlap(1)~=0 && dlap(2)==0 && dlap(3)==0 && dlap(4)~=0
                p4 = sind(teta(j,i))*slap(3);
                y4 = [p4*v(4)*(eps-dlap(4)) 0 p4*v(4)*dlap(4) -1 p4*v(4)];         % DARI 3 (Iso) ke 4 (Aniso)
                k4 = abs(asind(roots(y4)));
                if length(k4) > 1
                    for jj = 1:length(k4)
                        k44 = k4(jj)-teta(j,i);
                    end
                        [val,index] = min(k44);k4 = k4(index);
                end
                TT(j,i) = ((bounds(2)-bounds(1))/R(j))*RR(j,i)*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1)+...
                    ((bounds(3)-bounds(2))/R(j))*RR(j,i)*slap(3) + ((R(j)-bounds(3))/R(j))*RR(j,i)*(1+(-1*sind(k4)^2*cosd(k4)^2*dlap(4)))*slap(4);    
            elseif dlap(1)==0 && dlap(2)~=0 && dlap(3)~=0 && dlap(4)==0
                TT(j,i) = ((bounds(2)-bounds(1))/R(j))*RR(j,i)*(1+(-1*sind(k2)^2*cosd(k2)^2*dlap(2)))*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1)+...
                    ((bounds(3)-bounds(2))/R(j))*RR(j,i)*(1+(-1*sind(k23)^2*cosd(k23)^2*dlap(3)))*slap(3) + ((R(j)-bounds(3))/R(j))*RR(j,i)*slap(4);    
            elseif dlap(1)==0 && dlap(2)==0 && dlap(3)~=0 && dlap(4)~=0 || dlap(1)~=0 && dlap(2)==0 && dlap(3)~=0 && dlap(4)~=0
                p34 = sind(k3)*slap(3)*(1+(-1*sind(k3)^2*cosd(k3)^2*dlap(3)));
                y34 = [p34*v(4)*(eps-dlap(4)) 0 p34*v(4)*dlap(4) -1 p34*v(4)];             % DARI 2(Iso) ke 3 (Aniso) ke 4 (Aniso)
                k34 = abs(asind(roots(y34)));
                if length(k34) > 1
                    for jj = 1:length(k3)
                        k344 = k34(jj)-k3;
                    end
                        [val,index] = min(k344);k34 = k34(index);
                end
                TT(j,i) = ((bounds(2)-bounds(1))/R(j))*RR(j,i)*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1)+...
                    ((bounds(3)-bounds(2))/R(j))*RR(j,i)*(1+(-1*sind(k3)^2*cosd(k3)^2*dlap(3)))*slap(3) + ((R(j)-bounds(3))/R(j))*RR(j,i)*(1+(-1*sind(k34)^2*cosd(k34)^2*dlap(4)))*slap(4);  
            elseif dlap(1)==0 && dlap(2)~=0 && dlap(3)==0 && dlap(4)~=0 || dlap(1)~=0 && dlap(2)~=0 && dlap(3)==0 && dlap(4)~=0
                k230 =  asind(sind(k2)*slap(2)*(1+(-1*sind(k2)^2*cosd(k2)^2*dlap(2)))/slap(3));  %Teta Lapisan 3 (Iso), dari 2 (Aniso) ke 3 (Iso)
                p340 = sind(k230)*slap(3);
                y340 = [p340*v(4)*(eps-dlap(4)) 0 p340*v(4)*dlap(4) -1 p340*v(4)];         % DARI 2(Aniso) ke 3 (Iso) ke 4 (Aniso)
                k340 = abs(asind(roots(y340)));
                if length(k340) > 1
                    for jj = 1:length(k340)
                        k3341 = k340(jj)-k230;
                    end
                        [val,index] = min(k3341);k341 = k341(index);
                end
                TT(j,i) = ((bounds(2)-bounds(1))/R(j))*RR(j,i)*(1+(-1*sind(k2)^2*cosd(k2)^2*dlap(2)))*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1)+...
                    ((bounds(3)-bounds(2))/R(j))*RR(j,i)*slap(3) + ((R(j)-bounds(3))/R(j))*RR(j,i)*(1+(-1*sind(k340)^2*cosd(k340)^2*dlap(4)))*slap(4);
            elseif dlap(1)~=0 && dlap(2)~=0 && dlap(3)~=0 && dlap(4)==0
                TT(j,i) = ((bounds(2)-bounds(1))/R(j))*RR(j,i)*(1+(-1*sind(k2)^2*cosd(k2)^2*dlap(2)))*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1)+...
                    ((bounds(3)-bounds(2))/R(j))*RR(j,i)*(1+(-1*sind(k23)^2*cosd(k23)^2*dlap(3)))*slap(3) + ((R(j)-bounds(3))/R(j))*RR(j,i)*slap(4);
            elseif dlap(1)==0 && dlap(2)~=0 && dlap(3)~=0 && dlap(4)~=0 || dlap(1)~=0 && dlap(2)~=0 && dlap(3)~=0 && dlap(4)~=0
                p341 = sind(k23)*slap(3)*(1+(-1*sind(k23)^2*cosd(k23)^2*dlap(3)));
                y341 = [p341*v(4)*(eps-dlap(4)) 0 p341*v(4)*dlap(4) -1 p341*v(4)];         % DARI 2(Aniso) ke 3 (Aniso) ke 4 (Aniso)
                k341 = abs(asind(roots(y341)));
                if length(k341) > 1
                    for jj = 1:length(k341)
                        k3341 = k341(jj)-k23;
                    end
                        [val,index] = min(k3341);k341 = k341(index);
                end
                TT(j,i) = ((bounds(2)-bounds(1))/R(j))*RR(j,i)*(1+(-1*sind(k2)^2*cosd(k2)^2*dlap(2)))*slap(2) + (bounds(1)/R(j))*RR(j,i)*slap(1)+...
                    ((bounds(3)-bounds(2))/R(j))*RR(j,i)*(1+(-1*sind(k23)^2*cosd(k23)^2*dlap(3)))*slap(3) + ((R(j)-bounds(3))/R(j))*RR(j,i)*(1+(-1*sind(k341)^2*cosd(k341)^2*dlap(4)))*slap(4);
            end
            
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
    plot(shot(p),0,'O','MarkerSize',5,'MarkerFaceColor','r');hold on   
end
for p = 1:length(R)
    plot(0,R(p),'V','MarkerSize',8,'MarkerFaceColor','r');hold on     
end

for i=1:length(shot)
    for j =1:length(R)
        plot([shot(i),0],[0,R(j)],'-k');hold on     %PLOTTING RAY TRACING
    end
end
title('FORWARD MODELLING DOWNHOLE')
ylabel('DEPTH(m)');xlabel('OFFSET(m)')
grid on;set(gca,'ydir','reverse')