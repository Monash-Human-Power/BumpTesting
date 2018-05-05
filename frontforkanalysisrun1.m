clear all; close all; clc;

data = importdata('H3LIS331DL_front_run_1_2_26_pm.csv');


% Obtain time index
t = data.data(:,1);


%import acceleration data 

ax = 9.81*data.data(:,2); 
ay = 9.81*data.data(:,3); 
az = 9.81*data.data(:,4); 


% apply filter: Standards

fs = 1000;      % Sampling rate in Hz
T = 1/fs;
CFC = 166.62;
wd = 2 * pi * CFC * 2.0775;
wa = sin(wd * T/2) / cos(wd * T/2);
a0 = wa^2 / (1 + sqrt(2) * wa + wa^2);
a1 = 2 * a0;
a2 = a0;
b1 = -2 * (wa^2 - 1) / (1 + sqrt(2) * wa + wa^2);
b2 = (-1 + sqrt(2) * wa - wa^2) / (1 + sqrt(2) * wa + wa^2);

% Coefficient vectors
a = [a0 a1 a2];
b = [1 -b1 -b2];

% Display transfer function in command window
sys = tf(a,b);


% Filter Inputs
x = filtfilt(a,b,ax);
y = filtfilt(a,b,ay);
z = filtfilt(a,b,az);

maxaccelx = max(x)
maxaccely = max(y)
maxaccelz = max(z)

% plot raw data vs filterd data

figure
subplot(3,1,1)
plot(t,ax,'LineWidth',2)
hold on
plot(t,x,'LineWidth',1.2)
xlabel('time (s)')
ylabel('accel (m/s^2)')
legend('x-axis raw data','x-axis filtered')
subplot(3,1,2)
plot(t,ay,'LineWidth',2)
hold on
plot(t,y,'LineWidth',1.2)
xlabel('time (s)')
ylabel('accel (m/s^2)')
legend('y-axis raw data','y-axis filtered')
subplot(3,1,3)
plot(t,az,'LineWidth',2)
hold on
plot(t,z,'LineWidth',1.2)
xlabel('time (s)')
ylabel('accel (m/s^2)')
legend('z-axis raw data','z-axis filtered')
suptitle('raw and filtered acceleration vs time graphs')


%% integration
step = 0.0000001;
%find bounds manually with fitted function
% [fitx, pox] = PowerFitSplineDefault(t, x,0)
% curve fitted the data 'run1x' is the function
[fitx,pox] = run1x(t,x);


% bound where max oscilation occurs after peak but before 2nd zero
lbx = 3.604;
ubx = 3.61;

pt = [lbx:0.00001:ubx];
px = feval(fitx,pt);

plot(pt,px)

%if no zero choose min
% minqx = [4.293:0.00001:4.295];
% [mina,indexx] = min(feval(fitx,minqx))
% 
% qx = minqx(indexx)


 qx = fzero(fitx,[lbx,ubx])


% finding max acceleration between zeros
boundx = [qx:0.00001:ubx];
aboundx = feval(fitx,boundx);
[maxax,indx] = max(abs(aboundx))
% time where max a occurs
sx = boundx(indx);
%step for cumtrapz

% v and d integration from when a =0 to max a
vx = cumtrapz([qx:step:sx],[feval(fitx,[qx:step:sx])]);
dx = cumtrapz([qx:step:sx],vx);

if ubx == sx
    error('upper bound not right')
end


%%
[fity,poy] = run1y(t,y);
% bound where max oscilation occurs after peak but before 2nd zero
lby = 3.607;
uby = 3.61;

pty = [lby:0.00001:uby];
py = feval(fitx,pty);

figure
plot(pty,py)
%if no zero choose min
% minqy = [4.293:0.00001:4.295];
% [minay,indexy] = min(feval(fity,minqy))
% 
% qy = minqy(indexy)


qy = fzero(fity,[lby,uby])


% finding max acceleration between zeros
boundy = [qy:0.00001:uby];
aboundy = feval(fity,boundy);
[maxay,indy] = max(abs(aboundy))
% time where max a occurs
sy = boundy(indy);
%step for cumtrap

% v and d integration from when a =0 to max a
vy = cumtrapz([qy:step:sy],[feval(fity,[qy:step:sy])]);
dy = cumtrapz([qy:step:sy],vy);

if uby == sy
    error('upper bound not right')
end

%%
[fitz,poz] = run1z(t,z);
% bound where max oscilation occurs after peak but before 2nd zero
lbz = 3.428;
ubz = 3.431;

%if no zero choose min
% minqz = [4.293:0.00001:4.295];
% [minaz,indexz] = min(feval(fitz,minqz))
% 
% qz = minqz(indexz)


qz = fzero(fitz,[lbz,ubz])


% finding max acceleration between zeros
boundz = [qz:0.00001:ubz];
aboundz = feval(fitz,boundz);
[maxaz,indz] = max(abs(aboundz))
% time where max a occurs
sz = boundz(indz);
%step for cumtrapz


% v and d integration from when a =0 to max a
vz = cumtrapz([qz:step:sz],[feval(fitz,[qz:step:sz])]);
dz = cumtrapz([qz:step:sz],vz);

if ubz == sz
    error('upper bound not right')
end

%%

%plot velocity and displacement curves
figure
subplot(6,1,1)
plot([qx:step:sx],vx)
xlabel('time (s)')
ylabel('velocity (m/s)')
title('velocity x-axis')

subplot(6,1,2)
plot([qx:step:sx],dx)
xlabel('time (s)')
ylabel('displacement (m)')
title('displacement x-axis')

subplot(6,1,3)
plot([qy:step:sy],vy)
xlabel('time (s)')
ylabel('velocity (m/s)')
title('velocity y-axis')

subplot(6,1,4)
plot([qy:step:sy],dy)
xlabel('time (s)')
ylabel('displacement (m)')
title('displacement y-axis')

subplot(6,1,5)
plot([qz:step:sz],vz)
xlabel('time (s)')
ylabel('velocity (m/s)')
title('velocity z-axis')

subplot(6,1,6)
plot([qz:step:sz],dz)
xlabel('time (s)')
ylabel('displacement (m)')
title('displacement z-axis')

format short
qx
qy
qz
sx
sy
sz

format long
ex = dx(end)

ey = dy(end)

ez = dz(end)


