%TOMOGRAPHY SEISMIC
%WEAK ANISOTROPY
%KELOMPOK 1

clc;clear;close; 
%HATI-HATI KETIKA MEMASUKKAN v1,v2,v3,dan sudut awal KARENA MUDAH SEKALI LAPISAN 2 DAN 3 MENCAPAI SUDUT KRITIS 
v_lap = [750 900 1100] ; %kecepatan lapisan dalam m/s
h1 = 20;h2=60;h3=40; %Ketebalan lapisan
m  = linspace(2,40,10);           %Sudut datang lapisan pertama
p  = sind(m)./v_lap(1);  %Ray parameter
tl = [h1 h1+h2 h1+h2+h3]; %Lapisan Model
d  = 0.199;     %Parameter Anisotrop Faktor Sudut Fase
e  = 0.246;     %Fraksi perbedaan Vp vertikal terhadap horizontal
%Sudut Lapisan Kedua, NANTI masukkin persamaan ANISTROP  
%Persamaan anisotrop
%p*v2*(e-d)*(sind(teta))^4+p*v2*d*(sind(teta))^2-sind(teta)+p*v2 = 0
%tand(gamma)=tand(teta)*(1+2*d+4*(e-d)*(sind(teta)))   [GRUP VELO]

%MATRIX WADAH HASIL LOOPING
BIAS_2 = zeros(1,length(m));     %MATRIX SUDUT BIAS LAPISAN 2 PADA SUDUT DATANG n
BIAS_3 = BIAS_2;                 %MATRIX SUDUT BIAS LAPISAN 3 PADA SUDUT DATANG n
count = 0;

for i=1:length(m)
    a = p(i)*v_lap(2)*(e-d);l = p(i)*v_lap(2)*d;o = p(i)*v_lap(2);        %(111)
    y = @(x)a*x^4+l*x^2-x+o;
    yl = fzero(y,0.5);      %RUMUS NON-LINEAR MENDUGA-DUGA, KASIH TEBAKAN AWAL DULU (kali ini pakai 0.5)
    h = yl;
    krit =num2str(m(i));
    % y = [a 0 l -1 o];  %RUMUS SUDUT FASE
    % h = roots(y);
    H = asind(h);
    Y = atand(tand(H)*(1+2*d+4*(e-d)*(sind(H))^2)); %SUDUT GRUP LAPISAN ANISO
    YY = Y(isreal(Y)&Y<90&Y>0);
    BIAS_2(i) = YY;   %Check hasil real dan limit sudut 0<x<90 (kritis)
    m2 = BIAS_2(i);        %SUDUT BIAS LAPISAN 2 FIX
    % if isempty(m2)
    %     error('Tidak ada sudut yang real. Coba teliti di rumus fasenya')
    % end                                        %(222) 

    %COMMENT DARI (111) hingga (222)dan uncomment (333) untuk beralih ke model isotrop
    % m2 = asind(p*v2);               %(333)
    %COMMENT (333) dan uncomment (111)-(222) UNTUK BERALIH KE MODEL ANISOTROP
    if isreal(m2) ~= true
        fprintf('\nSudut bias lapisan kedua telah\nmelebihi batas kritis pada sudut datang %f \n',krit)
        break;
    elseif m2 == 90
        fprintf('\nSudut bias lapisan kedua telah\nmencapai kritis pada sudut datang %f \n',krit)
    end
    %Sudut Lapisan Ketiga\
    l = sind(m2)*v_lap(3)/v_lap(2);                   %SINUS SUDUT BIAS LAPISAN 3
    m3 = asind(sind(m2)*v_lap(3)/v_lap(2));           %SUDUT BIAS LAPISAN 3
    BIAS_3(i) = m3;
    if isreal(m2) == true && isreal(m3) ~= true 
        fprintf('\nSudut bias lapisan ketiga melebihi\nsudut kritis pada sudut datang %f derajat\n',krit)
    elseif m3 == 90
        fprintf('\nSudut bias lapisan ketiga telah\nmencapai kritis pada sudut datang %f derajat\n',krit)
    end
    %Lintasan rambat
    d1 = h1/cosd(m(i));
    d2 = h2/cosd(m2);
    d3 = h3/cosd(m3);
    %Time travel lintasan
    t1 = d1/v_lap(1);
    t2 = d2/v_lap(2);
    t3 = d3/v_lap(3);
    %TWT 
    T1 = 2*t1;
    T2 = T1+2*t2;
    T3 = T1+T2+2*t3;
    %Offset PermukaanT
    x1 = 2*tand(m(i))*h1;
    x2 = x1+2*tand(m2)*h2;
    x3 = x1+x2+2*tand(m3)*h3;
    %Titik Offset Refleksi
    X1 = [0.5*x1,h1]; %Reflek Border 1
    X2 = [0.5*x2,h2+h1]; %Refleks Border 2
    X3 = [0.5*x3,h3+h2+h1]; %Reflek Border 3
    X21= X1+[(2*tand(m2)*h2),0]; %Reflek Border 2 di Border 1
    X31= X1+[(2*(tand(m2)*h2+tand(m3)*h3)),0]; %Reflek Border 3 di Border 1 
    X32= X2+[(2*tand(m3)*h3),0]; %Reflek Border 3 di Border 2
    
    figure(1)
    %Plotting Model Block
    Offset = linspace(0,x3+10,50);
    Depth  = linspace(0,tl(3),50);
    plot(Offset,Depth)
    SDT = num2str(m(i));
    title('MODEL ANISOTROP SUDUT DATANG  Derajat')
    ylabel('Depth(m)')
    xlabel('Offset(m)')
    xlim([0 x3+10])
    ylim([0 tl(3)])
    yline(tl(1),'-k');hold on
    yline(tl(2),'-k');hold on
    xp =[0 0 x3+10 x3+10];
    yp =[0 tl(3) tl(3) 0];
    plot(xp,yp)
    Lx1 =[0 0 x3+10 x3+10];
    Ly1 =[0 tl(1) tl(1) 0];
    Ly2 =[tl(1) tl(2) tl(2) tl(1)];
    Ly3 =[tl(2) tl(3) tl(3) tl(2)];
    plot(Lx1,Ly1);               %LAPISAN ISOTROP 1
    fill(Lx1,Ly1,'y');          
    hold on
    plot(Lx1,Ly2);               %LAPISAN ANISOTROP
    fill(Lx1,Ly2,'g');
    hold on
    plot(Lx1,Ly3);               %LAPISAN ISOTROP 2
    fill(Lx1,Ly3,'b');
    hold on
    %Plotting Ray Tracing
    Hor_Refrak = [0 X1(1) X2(1) X3(1) x3-X2(1) x3-X1(1) x3];
    Dep_Refrak = [0 tl(1) tl(2) tl(3) tl(2) tl(1) 0];
    ReflekX_Border1 = [X1(1) x1];
    ReflekY_Border1 = [tl(1) 0];
    ReflekX_Border2 = [X2(1) X21(1) x2];
    ReflekY_Border2 = [tl(2) tl(1) 0];
    plot(Hor_Refrak,Dep_Refrak,'-r');hold on
    plot(ReflekX_Border1,ReflekY_Border1,'-r');hold on
    plot(ReflekX_Border2,ReflekY_Border2,'-r');hold on
    grid on
    set(gca,'ydir','reverse')
    figure(2)
    %TIME TRAVEL TABLE, Lnx = Data Offset Reflek ke-n, Lnt = Data Waktu Reflek ke-n
    L1x = [0 0.5*x1 x1];
    L1t = [0 t1 T1];
    L2x = [0 (0.5*x2-0.5*x1) (X21(1)-0.5*x1) (x2-0.5*x1)];
    L2t = [t1 t1+t2 t1+2*t2 T2];
    L3x = [0 (0.5*x3-0.5*x2) (X32(1)-0.5*x2) (X31(1)-0.5*x2) (x3-0.5*x2)];
    L3t = [t1+t2 t1+t2+t3 t1+t2+2*t3 t1+2*(t2+t3) T3];
    F=plot(L1x,L1t,'b');hold on
    G=plot(L2x,L2t,'g');hold on
    H=plot(L3x,L3t,'r');
    hold on
    legend([F,G,H],["Reflek 1","Reflek 2","Reflek 3"]);
    grid on
    title('TRAVEL TIME TABLE')
    xlabel('Offset(m)')
    ylabel('Time(s)')
    set(gca,'ydir','reverse')
    count = count+1;
end

%     if count~=length(m)
%         fprintf('Input sudut datang bermasalah pada data ke-%n',m(i));
%     end





