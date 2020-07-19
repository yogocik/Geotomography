%GEOTOMOGRAFI VSP 
%KELOMPOK 1

%% PARAMETER AWAL
clear;clc;close;

teta1   = 60; %Sudut datang untuk semua shot dari lapisan pertama
L       = 6;  %Lapisan Tanah
v       = zeros(L,1); % Vektor V tiap lapisan;
a       = 400; % V lapisan pertama
v_inc   = 20; % Kenaikan kecepatan tiap lapisan
rec     = [40,100,170,250]; % Shot Offset
tbl     = 20; %tebal lapisan pertama
for i=1:L
    if i == 1
        v(i) = a;
    else
        v(i) = v(i-1) + v_inc;
    end
end


%RECEIVER BERADA DI TENGAH LAPISAN
s = sind(teta1)/v(1);      %SLOWNESS TIAP LAPISAN (dari lapisan pertama)
r = length(L) ;              %JUMLAH RECEIVER
sh= length(rec) ;              %JUMLAH SHOT 

%SUDUT SHOT 2 
d1  = s*v(2);
BL1 = asin(d1);
%SUDUT SHOT 3
d2  =s*v(3);
BL2 = asind(d2);
%SUDUT SHOT 4
d3 =s*v(4);
BL3 = asind(d3);

%JARAK VERTIKAL RECEIVER 1
JV1 = 0.5*tbl; % Posisi R1 dari titik 0 m
%JARAK VERTIKAL RECEIVER 2
x1 = tbl*tand(teta1); 
xx1 = rec(2)-x1; 
delta_H = tand(BL1)*xx1; %Posisi R2 dari BL1
JV2 = 2*JV1 + delta_H; %Posisi R2 dari titik 0 m
tbl2 = 2*delta_H; % Tebal Lapisan 2
%JARAK VERTIKAL RECEIVER 3
xxx1 =rec(3)-x1-tbl2*tand(BL2); 
delta_H2 = tand(BL2)*xxx1;%posisi R3 dari BL2
tbl3 = 2*delta_H2; %tebal lapisan 3
JV3 = 2*JV1 + tbl2 +delta_H2; %Posisi R3 dari titik 0 m
%JARAK VERTIKAL RECEIVER 4
xxxx1 = rec(4)-x1-tbl2*tand(BL2)-tbl3*tand(BL3);
delta_H3 = tand(BL3)*xxxx1; % Posisi R4 dari BL3
tbl4 = 2*delta_H3;   % Tebal lapisan 4
JV4 =2*JV1 + tbl2 + tbl3 + delta_H3; % Posisi R4 dari titik 0 m


        
        








      