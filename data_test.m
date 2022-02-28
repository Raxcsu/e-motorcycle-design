clc; clear all; close all;
P = gpxread('daniel_27_02_2022.gpx');

Lat = P.Latitude;
Long = P.Longitude;
T = P.Time;

mv = 110; %masa del vehiculo kg
mp = 90; %masa persona kg 
me = 10; %masa del compartimiento kg

M = mv + mp + me; % 15.2/19.2

E2 = 0; %energia 2

Time = [];

for i =1:length(Long)
    
   temp = deg2km(Lat(i))*10e+02;
   X(i) = temp;
   a = temp;
   temp = deg2km(Long(i))*10e+02;
   Y(i) = temp;
   b = temp; 
   
   vec(i) = sqrt(a^2+b^2);
   
   temp = T{i};
   horas = temp(13);
   horas = str2num(horas);
   
   minutos = temp(15);
   minutos = str2num(minutos)*10;
   
   minutos2 = temp(16);
   minutos2 = str2num(minutos2);
   
   minutos = minutos2 + minutos;
   
   sec = temp(18);
   sec = str2num(sec)*10;
   
   sec2 = temp(19);
   sec2 = str2num(sec2);
   
   sec = sec2 + sec;
   
   temp = horas*3600 + minutos*60 + sec;  %tiempo en segundos.
   Time(i) = temp;
end

vel = [];
velm = [];
velm(1) = 0;
vel(1) = 0;

d = 0; % distancia recorrida
Disp = [];
Disp(1) = 0;

accel = [];
accel(1) = 0;

for i=2:length(Long)
    
   	Dx = X(i) - X(i-1);
    Dy = Y(i) - Y(i-1);
    Dt = Time(i) - Time(i-1);
    temp = sqrt(Dx^2 + Dy^2);
    Disp(i) = temp;
    d = d + temp;
    temp = temp/Dt;
    velm(i) = temp;
    temp = temp*3.6; %convert to kmph
    vel(i) = temp;
    
    temp = velm(i) - velm(i-1);
    temp = temp/Dt;
    accel(i) = temp;
    
end 

%calculo de velocidad promedio y energias
velp = 0;
n = length(vel);
W = 0;
%accel = movmean(accel,3)

for i = 1:n %length(vel)
    
    if isnan(vel(i)) | (vel(i) == Inf)
         vel(i)
    else
        velp = velp + vel(i);
    end
    
    if ~(isnan(vec(i)) | (isnan(accel(i))) | abs(accel(i)) == Inf | abs(vec(i)) == Inf)
        if (accel(i) > 0 && accel(i) < 20)
             E2 = E2 + M*(accel(i))*Disp(i); %calculo de energia. 
        end 
    end 
end  

disp('Velocidad promedio (km)')
velp = velp/n %length(vel)
%velocidad promedio
disp('Energia (kWh)')
E2 = E2/3.6e+6 %convertimos de J a kWh 
disp('Distancia (km)')
d = d/1000 %distancia en km
subplot(3,1,1);
plot(Time,vec);grid;xlabel('Tiempo[s]'),ylabel('Desplazamiento [m]');

subplot(3,1,2); 
plot(Time,movmean(vel,30));grid;xlabel('Tiempo[s]'),ylabel('Velocidad kmph');

subplot(3,1,3); 
plot(Time,accel);grid;xlabel('Tiempo[s]'),ylabel('acceleracion m/s2');

%%
%calculation of the average velocity without considering standby times.
x = input('input range (km) desired: ');
velac = 0;
p = 0;
mv = 7; % minimum velocity
for i = 1:length(vel)
   if (vel(i) > mv) && ~(isnan(vel(i)) | (vel(i) == Inf))
       p = p + 1; 
       velac = velac + vel(i);
   end 
end
disp('average velocity (kmph)')
vpf = velac/p

% efficiency: km/kwh
disp('km /kwh ')
ef = d/E2
disp(' Capacity required (kWh): ')
C = x/ef

%%
m1 = 0;
for i = 1:length(accel)
    if accel(i) >= 0
        m1 = m1 + 1;
    end 
end 

r = m1/length(accel)
figure(2),pie(r)
