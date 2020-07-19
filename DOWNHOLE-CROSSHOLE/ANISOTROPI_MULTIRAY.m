%% TOMOGRAPHY SEISMIC REFLECTION REV 2.1
%FORWARD MODELLING : SIMPLE WEAK ANISOTROPY 3 LAYERS (ONE SOURCE, MULTI-RAYPATH)
%KELOMPOK 1
%ETS 2019/2020
%<--------SHORT INFO--------->
%MOHON PERGUNAKAN SCRIPT SEBIJAK-BIJAKNYA DENGAN TIDAK MENGUBAH ALGORITMA
%KECUALI NILAI INPUT BERUPA
%v_lap,h,angle_lower,angle_upper,jumlah_rays,d,e,L dan T.
%SEMOGA SCRIPT BISA BERMANFAAT UNTUK BELAJAR RAY TRACING DAN TRAVEL TIMES
%MODEL ANISOTROPY LEMAH

%INSTITUT TEKNOLOGI SEPULUH NOPEMBER----23 OKTOBER 2019
%DEPARTEMEN TEKNIK GEOFISIKA
%FAKULTAS TEKNIK SIPIL,LINGKUNGAN, DAN KEBUMIAN

clc;clear;close

%HATI-HATI KETIKA MEMASUKKAN v1,v2,v3,dan sudut awal... 
%KARENA MUDAH SEKALI LAPISAN 2 DAN 3 MENCAPAI SUDUT KRITIS
%JANGAN MEMBUAT NILAI v_lap DENGAN PERBEDAAN YANG BESAR
v_lap    = [600 900 1200] ;     %Kecepatan lapisan dalam m/s
h        = [10 20 30];          %Ketebalan lapisan
angle_lower = 5;   %BATAS BAWAH SUDUT DATANG AWAL         
angle_upper = 25;  %BATAS ATAS SUDUT DATANG AWAL
jumlah_rays = 60;   %JUMLAH RAY YANG DIINGINKAN
m        = linspace(angle_lower,angle_upper,jumlah_rays);   %Sudut datang lapisan pertama
p        = sind(m)./v_lap(1);                               %Ray parameter
tl       = [h(1) h(1)+h(2) h(1)+h(2)+h(3)];                 %Lapisan Model
d        = 0.199;               %Parameter Anisotrop Faktor Sudut Fase
e        = 0.246;               %Fraksi perbedaan Vp vertikal terhadap horizontal                           
L =0.1;                         %ROOTS SOLVER INITIATOR (fzero parameter)
T =2;                           %Posisi Shot/SOURCE
%MATRIX WADAH HASIL LOOPING
PHASE  = zeros(size(m));         %SUDUT FASE LAPISAN ANISOTROPI
GROUP  = PHASE;                  %SUDUT GROUP LAPISAN ANISTROPI
count  = 0;                      %ITERASI PROGRESS 
k = 0;q = 0;s = 0;               %SHELL

%PERSAMAAN ANISOTROPI (TIDAK USAH DIBACA DARIPADA BINGUNG!!)
%p*v2*(e-d)*(sind(teta))^4+p*v2*d*(sind(teta))^2-sind(teta)+p*v2=0 [PHASE]
%tand(gamma)=tand(teta)*(1+2*d+4*(e-d)*(sind(teta)))   [GROUP]
%LAPISAN KEDUA BERSIFAT ANISOTROPI
fprintf('\n1.KETEBALAN LAPISAN 1,2,3 YAITU %d METER,%d METER,dan %d METER\n',h)
fprintf('\n2.BATAS LAPISAN 1 DAN 2 BERADA PADA KEDALAMAN %d METER dan %d METER\n',[tl(1) tl(2)]) 
%PERHITUNGAN SUDUT LAPISsAN ANISOTROPI
for i=1:length(m)
    y = @(x) p(i)*v_lap(2)*(e-d)*(x)^4+p(i)*v_lap(2)*d*(x)^2-x+p(i)*v_lap(2);
    k = asind(fzero(y,L));
    if isreal(k)~= true
        fprintf('\n3.1.SUDUT FASE LAPISAN KEDUA MELEBIHI BATAS KRITIS\nPADA SUDUT DATANG: %f DERAJAT\n',m(i))
        break
    elseif k == 90
        fprintf('\n3.2.SUDUT FASE LAPISAN KEDUA MENCAPAI KRITIS\nPADA SUDUT DATANG: %f DERAJAT',m(i))
    end
    PHASE(i) = k;
    q        = atand(tand(k)*(p(i)*v_lap(2)*d+2*d+4*(e-d)*(sind(k))));
    if isreal(q)~= true
        fprintf('\n4.1.SUDUT GRUP LAPISAN KEDUA MELEBIHI BATAS KRITIS\n PADA SUDUT AWAL: %f DERAJAT\n',m(i))
        break
    elseif q == 90
        fprinft('\4.1.1nSUDUT GRUP LAPISAN KEDUA MENCAPAI KRITIS\n PADA SUDUT AWAL: %f DERAJAT\n',m(i))
    end
    GROUP(i) = q;
    count = count+1;
end

GROUP = GROUP(GROUP~=0)';      %SUDUT GRUP FIX
PHASE = PHASE(PHASE~=0)';      %SUDUT FASE FIX
BIAS_3 = zeros(size(GROUP));   %SHELL SUDUT LAPISAN KETIGA
COUNT = 0;

%PERHITUNGAN SUDUT LAPISAN KETIGA (ISOTROP)
for j =1:length(BIAS_3)
    s = asind(sind(GROUP(j))*v_lap(3)/v_lap(2));
    if isreal(s)~=true
        fprintf('\n5.1.SUDUT BIAS LAPISAN KETIGA MELEBIHI BATAS KRITIS\n PADA SUDUT BIAS KEDUA: %f DERAJAT\n',m(j))
        break;
    elseif s == 90
        fprintf('\n5.2.SUDUT BIAS LAPISAN KETIGA MENCAPAI BATAS KRITIS\n PADA SUDUT BIAS KEDUA: %f DERAJAT\n',m(j))
    end
   BIAS_3(j) = s;
   COUNT = COUNT+1;
end

BIAS_3 = BIAS_3(BIAS_3~=0);    %MATRIX SUDUT LAPISAN KETIGA FIX
fprintf('\n6.PERHITUNGAN LAPISAN ANISOTROPI\nSAMPAI PADA INPUT DATA SUDUT AWAL KE-%d\n',count)
fprintf('\n7.PERHITUNGAN LAPISAN KETIGA\nSAMPAI PADA INPUT DATA SUDUT BIAS KEDUA KE-%d\n',COUNT)
fprintf('\n8.RAY TRACING FULL DIHASILKAN\nHINGGA DATA SUDUT AWAL KE-%d\n',COUNT) 

%SORTING DATA SUDUT
BIAS_2_FIX =  GROUP(1:length(BIAS_3)); %MATRIX SUDUT LAPISAN ANISOTROP YANG DIPAKAI
SUDUT = [m(1:length(BIAS_2_FIX))' BIAS_2_FIX BIAS_3]; %MATRIX SEMUA SUDUT 
RAY = zeros(size(SUDUT));                             %MATRIX SHELL RAYPATH
LON = RAY;                                            %MATRIX SHELL OFFSET PER LAPISAN
for j = 1:size(SUDUT,1)                               %PERUBAHAN POLA RAY TRACING
    if v_lap(2)>v_lap(1)
        if SUDUT(j,2)<SUDUT(j,1)
            fprintf('\n9.a.POLA RAY TRACING LAPISAN ANISOTROPI MENDEKATI GARIS NORMAL\nPADA DATA SUDUT AWAL: %f DERAJAT (v_aniso(lapisan_2)>v_lapisan_11)\n',m(j))
            fprintf('\n9.b.PERUBAHAN POLA RAY TRACING DIPENGARUHI OLEH INPUT PARAMETER fzero,\nINTERVAL SUDUT AWAL,DAN JUMLAH RAYS\n') 
            break;
        end
    end    
end

%PERHITUNGAN RAYPATH
for j = 1:length(v_lap) %PER LAPISAN
    for i=1:length(BIAS_2_FIX) %Per Ray/Sudut Awal
        RAY(i,j) = tand(SUDUT(i,j))*h(j);        %PANJANG OFFSET PER LAPIS SEMUA SUDUT
        LON(i,j) = (1/cosd(SUDUT(i,j)))*h(j);    %PANJANG RAYPATH PER LAPIS PER RAY
    end
end

RAY_TOTAL   = zeros(size(RAY));       %KOORDINAT TITIK AKUMULASI OFFSET REFLEKSI TIAP BATAS LAPIS 
Travel_Time = zeros(size(RAY_TOTAL)); %WAKTU TEMPUH PER LAPIS PER SUDUT

for j = 1:length(BIAS_2_FIX)
    RAY_TOTAL(j,1) = RAY(j,1)+T;                %OFFSET AKUMULATIF (SEBENARNYA)
    RAY_TOTAL(j,2) = RAY(j,1) + RAY(j,2)+T;
    RAY_TOTAL(j,3) = RAY(j,1) + RAY(j,2) + RAY(j,3)+T;
    Travel_Time(j,1) = LON(j,1)/v_lap(1);
    Travel_Time(j,2) = LON(j,2)/v_lap(2);
    Travel_Time(j,3) = LON(j,3)/v_lap(3);
end

TRAV_LPN_1 = Travel_Time(:,1);  %LAPISAN 1 ALL SUDUT AWAL
TRAV_LPN_2 = Travel_Time(:,2);  %LAPISAN 2 ALL SUDUT AWAL
TRAV_LPN_3 = Travel_Time(:,3);  %LAPISAN 3 ALL SUDUT AWAL


%TITIK REFLEKSI PER LAPISAN
LPN_1   = [RAY_TOTAL(:,1) ones(size(BIAS_2_FIX)).*tl(1)];
LPN_2   = [RAY_TOTAL(:,2) ones(size(BIAS_2_FIX)).*tl(2)];
LPN_3   = [RAY_TOTAL(:,3) ones(size(BIAS_2_FIX)).*tl(3)];
LPN_10  = [(RAY_TOTAL(:,1)+RAY(:,1)) zeros(size(BIAS_2_FIX))];
LPN_21  = [(RAY_TOTAL(:,2)+RAY(:,2)) ones(size(BIAS_2_FIX)).*tl(1)];
LPN_20  = [(RAY_TOTAL(:,2)+RAY(:,2)+RAY(:,1)) zeros(size(BIAS_2_FIX))];
LPN_32  = [(RAY_TOTAL(:,3)+RAY(:,3)) ones(size(BIAS_2_FIX)).*tl(2)];
LPN_31  = [(RAY_TOTAL(:,3)+RAY(:,3)+RAY(:,2)) ones(size(BIAS_2_FIX)).*tl(1)];
LPN_30  = [(RAY_TOTAL(:,3)+RAY(:,3)+RAY(:,2)+RAY(:,1)) zeros(size(BIAS_2_FIX))];

fprintf('\n10.JIKA BERMASALAH,COMMENT LINE 169-170 (axis limit)\n')
fprintf('\n11.JIKA MASIH BERMASALAH,COMMENT LINE 178 (legend function)\n')

%PLOTTING

figure(2)
X  = [ones(size(BIAS_2_FIX)).*T zeros(size(BIAS_2_FIX))];      %MATRIX POSISI SOURCE
XX = max(RAY_TOTAL);
TT = (ones(size(RAY)).*tl);                                    %MATRIX BATAS LAPISAN
plot(T,zeros(length(T)),'xr','MarkerSize',12);hold on          %TITIK SHOT
for i=1:length(tl)
    plot(linspace(0,2*XX(3),length(BIAS_2_FIX)),TT(:,i),'--k');hold on
end
for p = 1:length(BIAS_2_FIX)                                   %RAY TRACING
    plot([X(p,1) LPN_1(p,1)],[X(p,2) LPN_1(p,2)]);hold on
    plot([LPN_1(p,1) LPN_2(p,1)],[LPN_1(p,2) LPN_2(p,2)]);hold on
    plot([LPN_2(p,1) LPN_3(p,1)],[LPN_2(p,2) LPN_3(p,2)]);hold on
    plot([LPN_3(p,1) LPN_32(p,1)],[LPN_3(p,2) LPN_32(p,2)]);hold on
    plot([LPN_32(p,1) LPN_31(p,1)],[LPN_32(p,2) LPN_31(p,2)]);hold on
    plot([LPN_31(p,1) LPN_30(p,1)],[LPN_31(p,2) LPN_30(p,2)]);hold on
    plot([LPN_1(p,1) LPN_10(p,1)],[LPN_1(p,2) LPN_10(p,2)]);hold on
    plot([LPN_2(p,1) LPN_21(p,1)],[LPN_2(p,2) LPN_21(p,2)]);hold on
    plot([LPN_21(p,1) LPN_20(p,1)],[LPN_21(p,2) LPN_20(p,2)]);hold on
end
title('RAY TRACING REFLEKSI')
xlabel('Offset(m)')
ylabel('Depth(m)')
ylim([0 tl(3)]);hold on
xlim([0 2*XX(3)])
grid on
set(gca,'ydir','reverse')

figure(3)
R1=plot(RAY(:,1),TRAV_LPN_1,'-b');hold on
R2=plot(RAY(:,2),TRAV_LPN_2,'-g');hold on
R3=plot(RAY(:,3),TRAV_LPN_3,'-r');hold on
legend([R1,R2,R3],["Refleksi 1","Refleksi 2","Refleksi 3"])
grid on
title('TRAVEL TIME REFLEKSI')
ylabel('Waktu RayPath(s)')
xlabel('Offset RayPath(m)')
xlim([eps*-1 max(RAY(:,3))+5])
set(gca,'ydir','reverse')






    