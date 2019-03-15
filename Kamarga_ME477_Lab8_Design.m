% Khrisna Kamarga
% ME477 - Lab7
clear all; close all; clc;

k_vi = 0.41;
k_t = 0.11;
J = 3.8e-4;
B = 0;
T = 5e-3;
k = k_vi*k_t;

plant = zpk([], [-B/J 0], k/J);
s = tf('s');
fb = 8; % control bandwidth (Hz)
wc = 2*pi*fb;

options = pidtuneOptions('DesignFocus', 'reference-tracking');
c = pidtune(plant, 'pidf', wc, options);

G = feedback(series(c,plant), 1);
H = k*G/(plant); %torque
step(G)
figure(2)
step(H)

cdp = c2d(c, T, 'Tustin');
cd = tf(cdp);
[b a] = tfdata(cd, 'v');
sos = tf2sos(b, a);

HeaderFileName = 'C:\Users\Khrisna Adi Kamarga\Desktop\ME477\myPIDF.h';
fid = fopen(HeaderFileName, 'W');
comment = 'DC Motor Position Controller';
sos2header(fid, sos, 'PIDF', T, comment);
fclose(fid);


%%
load Lab8_Khrisna
close all; clc;
T = 5e-3;
t = T*(1:length(pref));
Gd = feedback(series(cd, c2d(plant,T)), 1);
Hd = Gd/(c2d(plant, T));
y = lsim(Gd, pref, t);
y_torq = 2*pi*k*lsim(Hd, pref, t);


subplot(3,1,1)
hold on
plot(t, y)
plot(t, pact, '.')
plot(t, pref)
title("Theoretical, Actual, and Reference Position");
legend theoretical actual reference
hold off

subplot(3,1,2)
hold on
plot(t, y - pref)
plot(t, pact - pref)
title("Theoretical and Actual Position Error")
legend theoretical actual
hold off

subplot(3,1,3)
hold on
plot(t, torq/1000)
plot(t, y_torq)
title("Theoretical and Actual Torque")
legend actual theoretical


