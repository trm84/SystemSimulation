%Tyler Matthews 4/15/19
%System Simulation Problem 13

close all; clc; clear; %Clear Console and Close Figures

%% Alpha = 0.1 - Initial Condition x,y,z=0.1

N = 10000;
T = 0.1;
t = linspace(0, T*N, N+1);

alpha = 0.1;

fxc = zeros(1,N);
fyc = zeros(1,N);
fzc = zeros(1,N);

xp = zeros(1,N);
yp = zeros(1,N);
zp = zeros(1,N);

fxp = zeros(1,N);
fyp = zeros(1,N);
fzp = zeros(1,N);

xc = zeros(1,N);
yc = zeros(1,N);
zc = zeros(1,N);

% Initial Conditions %
xc(1) = 0.1;
yc(1) = 0.1;
zc(1) = 0.1;

for k = 1:N
    % Evaluate for prediction %
    fxc(k) = alpha*(yc(k)+((xc(k)-2*(xc(k)).^3)/(7)));
    fyc(k) = xc(k)-yc(k)+zc(k);
    fzc(k) = -(100/7)*yc(k);
    % Predict  - Euler's Method
    xp(k+1) = xc(k) + T*fxc(k); %x(n+1) = x(n) + T*x'(n)
    yp(k+1) = yc(k) + T*fyc(k);
    zp(k+1) = zc(k) + T*fzc(k);
    
    % Evaluate for correction
    fxp(k+1) = alpha*(yp(k+1)+((xp(k+1)-2*(xp(k+1)).^3)/(7)));
    fyp(k+1) = xp(k+1)-yp(k+1)+zp(k+1);
    fzp(k+1) = -(100/7)*yp(k+1);
    % Correct Trapezoidal Method
    xc(k+1) = xc(k) + (T/2)*(fxp(k+1) + fxp(k)); %x(n+1) = x(n) + (T/2) * (x'(n+1) +x'(n))
    yc(k+1) = yc(k) + (T/2)*(fyp(k+1) + fyp(k));
    zc(k+1) = zc(k) + (T/2)*(fzp(k+1) + fzp(k));
end

figure;
subplot(2,2,1)
plot3(xc,yc,zc)
title('Chuas System: 3D Plot')


subplot(2,2,2)
plot(xc, t)
ylim([0,10])
title('X vs T')

subplot(2,2,3)
plot(yc, t)
ylim([0,10])
title('Y vs T')

subplot(2,2,4)
plot(zc, t)
ylim([0,10])
title('Z vs T')








%% Alpha = 1 - Initial Condition x,y,z=0.1
clear; %clear workspace

N = 10000;
T = 0.1;
t = linspace(0, T*N, N+1);

alpha = 1;

fxc = zeros(1,N);
fyc = zeros(1,N);
fzc = zeros(1,N);

xp = zeros(1,N);
yp = zeros(1,N);
zp = zeros(1,N);

fxp = zeros(1,N);
fyp = zeros(1,N);
fzp = zeros(1,N);

xc = zeros(1,N);
yc = zeros(1,N);
zc = zeros(1,N);

% Initial Conditions %
xc(1) = 0.1;
yc(1) = 0.1;
zc(1) = 0.1;

for k = 1:N
    % Evaluate for prediction %
    fxc(k) = alpha*(yc(k)+((xc(k)-2*(xc(k)).^3)/(7)));
    fyc(k) = xc(k)-yc(k)+zc(k);
    fzc(k) = -(100/7)*yc(k);
    % Predict  - Euler's Method
    xp(k+1) = xc(k) + T*fxc(k); %x(n+1) = x(n) + T*x'(n)
    yp(k+1) = yc(k) + T*fyc(k);
    zp(k+1) = zc(k) + T*fzc(k);
    
    % Evaluate for correction
    fxp(k+1) = alpha*(yp(k+1)+((xp(k+1)-2*(xp(k+1)).^3)/(7)));
    fyp(k+1) = xp(k+1)-yp(k+1)+zp(k+1);
    fzp(k+1) = -(100/7)*yp(k+1);
    % Correct Trapezoidal Method
    xc(k+1) = xc(k) + (T/2)*(fxp(k+1) + fxp(k)); %x(n+1) = x(n) + (T/2) * (x'(n+1) +x'(n))
    yc(k+1) = yc(k) + (T/2)*(fyp(k+1) + fyp(k));
    zc(k+1) = zc(k) + (T/2)*(fzp(k+1) + fzp(k));
end

figure;
subplot(2,2,1)
plot3(xc,yc,zc)
title('Chuas System: 3D Plot')


subplot(2,2,2)
plot(xc, t)
ylim([0,10])
title('X vs T')

subplot(2,2,3)
plot(yc, t)
ylim([0,10])
title('Y vs T')

subplot(2,2,4)
plot(zc, t)
ylim([0,10])
title('Z vs T')




%% Alpha = 10 - Initial Condition x,y,z=0.1
clear; %clear workspace

N = 10000;
T = 0.1;
t = linspace(0, T*N, N+1);

alpha = 10;

fxc = zeros(1,N);
fyc = zeros(1,N);
fzc = zeros(1,N);

xp = zeros(1,N);
yp = zeros(1,N);
zp = zeros(1,N);

fxp = zeros(1,N);
fyp = zeros(1,N);
fzp = zeros(1,N);

xc = zeros(1,N);
yc = zeros(1,N);
zc = zeros(1,N);

% Initial Conditions %
xc(1) = 0.1;
yc(1) = 0.1;
zc(1) = 0.1;

for k = 1:N
    % Evaluate for prediction %
    fxc(k) = alpha*(yc(k)+((xc(k)-2*(xc(k)).^3)/(7)));
    fyc(k) = xc(k)-yc(k)+zc(k);
    fzc(k) = -(100/7)*yc(k);
    % Predict  - Euler's Method
    xp(k+1) = xc(k) + T*fxc(k); %x(n+1) = x(n) + T*x'(n)
    yp(k+1) = yc(k) + T*fyc(k);
    zp(k+1) = zc(k) + T*fzc(k);
    
    % Evaluate for correction
    fxp(k+1) = alpha*(yp(k+1)+((xp(k+1)-2*(xp(k+1)).^3)/(7)));
    fyp(k+1) = xp(k+1)-yp(k+1)+zp(k+1);
    fzp(k+1) = -(100/7)*yp(k+1);
    % Correct Trapezoidal Method
    xc(k+1) = xc(k) + (T/2)*(fxp(k+1) + fxp(k)); %x(n+1) = x(n) + (T/2) * (x'(n+1) +x'(n))
    yc(k+1) = yc(k) + (T/2)*(fyp(k+1) + fyp(k));
    zc(k+1) = zc(k) + (T/2)*(fzp(k+1) + fzp(k));
end

figure;
subplot(2,2,1)
plot3(xc,yc,zc)
title('Chuas System: 3D Plot')


subplot(2,2,2)
plot(xc, t)
ylim([0,10])
title('X vs T')

subplot(2,2,3)
plot(yc, t)
ylim([0,10])
title('Y vs T')

subplot(2,2,4)
plot(zc, t)
ylim([0,10])
title('Z vs T')

%% Alpha = 0.1 - Initial Condition x,y,z=1.5

N = 10000;
T = 0.1;
t = linspace(0, T*N, N+1);

alpha = 0.1;

fxc = zeros(1,N);
fyc = zeros(1,N);
fzc = zeros(1,N);

xp = zeros(1,N);
yp = zeros(1,N);
zp = zeros(1,N);

fxp = zeros(1,N);
fyp = zeros(1,N);
fzp = zeros(1,N);

xc = zeros(1,N);
yc = zeros(1,N);
zc = zeros(1,N);

% Initial Conditions %
xc(1) = 1.5;
yc(1) = 1.5;
zc(1) = 1.5;

for k = 1:N
    % Evaluate for prediction %
    fxc(k) = alpha*(yc(k)+((xc(k)-2*(xc(k)).^3)/(7)));
    fyc(k) = xc(k)-yc(k)+zc(k);
    fzc(k) = -(100/7)*yc(k);
    % Predict  - Euler's Method
    xp(k+1) = xc(k) + T*fxc(k); %x(n+1) = x(n) + T*x'(n)
    yp(k+1) = yc(k) + T*fyc(k);
    zp(k+1) = zc(k) + T*fzc(k);
    
    % Evaluate for correction
    fxp(k+1) = alpha*(yp(k+1)+((xp(k+1)-2*(xp(k+1)).^3)/(7)));
    fyp(k+1) = xp(k+1)-yp(k+1)+zp(k+1);
    fzp(k+1) = -(100/7)*yp(k+1);
    % Correct Trapezoidal Method
    xc(k+1) = xc(k) + (T/2)*(fxp(k+1) + fxp(k)); %x(n+1) = x(n) + (T/2) * (x'(n+1) +x'(n))
    yc(k+1) = yc(k) + (T/2)*(fyp(k+1) + fyp(k));
    zc(k+1) = zc(k) + (T/2)*(fzp(k+1) + fzp(k));
end

figure;
subplot(2,2,1)
plot3(xc,yc,zc)
title('Chuas System: 3D Plot')


subplot(2,2,2)
plot(xc, t)
ylim([0,10])
title('X vs T')

subplot(2,2,3)
plot(yc, t)
ylim([0,10])
title('Y vs T')

subplot(2,2,4)
plot(zc, t)
ylim([0,10])
title('Z vs T')








%% Alpha = 1 - Initial Condition x,y,z=1.5
clear; %clear workspace

N = 10000;
T = 0.1;
t = linspace(0, T*N, N+1);

alpha = 1;

fxc = zeros(1,N);
fyc = zeros(1,N);
fzc = zeros(1,N);

xp = zeros(1,N);
yp = zeros(1,N);
zp = zeros(1,N);

fxp = zeros(1,N);
fyp = zeros(1,N);
fzp = zeros(1,N);

xc = zeros(1,N);
yc = zeros(1,N);
zc = zeros(1,N);

% Initial Conditions %
xc(1) = 1.5;
yc(1) = 1.5;
zc(1) = 1.5;

for k = 1:N
    % Evaluate for prediction %
    fxc(k) = alpha*(yc(k)+((xc(k)-2*(xc(k)).^3)/(7)));
    fyc(k) = xc(k)-yc(k)+zc(k);
    fzc(k) = -(100/7)*yc(k);
    % Predict  - Euler's Method
    xp(k+1) = xc(k) + T*fxc(k); %x(n+1) = x(n) + T*x'(n)
    yp(k+1) = yc(k) + T*fyc(k);
    zp(k+1) = zc(k) + T*fzc(k);
    
    % Evaluate for correction
    fxp(k+1) = alpha*(yp(k+1)+((xp(k+1)-2*(xp(k+1)).^3)/(7)));
    fyp(k+1) = xp(k+1)-yp(k+1)+zp(k+1);
    fzp(k+1) = -(100/7)*yp(k+1);
    % Correct Trapezoidal Method
    xc(k+1) = xc(k) + (T/2)*(fxp(k+1) + fxp(k)); %x(n+1) = x(n) + (T/2) * (x'(n+1) +x'(n))
    yc(k+1) = yc(k) + (T/2)*(fyp(k+1) + fyp(k));
    zc(k+1) = zc(k) + (T/2)*(fzp(k+1) + fzp(k));
end

figure;
subplot(2,2,1)
plot3(xc,yc,zc)
title('Chuas System: 3D Plot')


subplot(2,2,2)
plot(xc, t)
ylim([0,10])
title('X vs T')

subplot(2,2,3)
plot(yc, t)
ylim([0,10])
title('Y vs T')

subplot(2,2,4)
plot(zc, t)
ylim([0,10])
title('Z vs T')




%% Alpha = 10 - Initial Condition x,y,z=1.5
clear; %clear workspace

N = 10000;
T = 0.1;
t = linspace(0, T*N, N+1);

alpha = 10;

fxc = zeros(1,N);
fyc = zeros(1,N);
fzc = zeros(1,N);

xp = zeros(1,N);
yp = zeros(1,N);
zp = zeros(1,N);

fxp = zeros(1,N);
fyp = zeros(1,N);
fzp = zeros(1,N);

xc = zeros(1,N);
yc = zeros(1,N);
zc = zeros(1,N);

% Initial Conditions %
xc(1) = 1.5;
yc(1) = 1.5;
zc(1) = 1.5;

for k = 1:N
    % Evaluate for prediction %
    fxc(k) = alpha*(yc(k)+((xc(k)-2*(xc(k)).^3)/(7)));
    fyc(k) = xc(k)-yc(k)+zc(k);
    fzc(k) = -(100/7)*yc(k);
    % Predict  - Euler's Method
    xp(k+1) = xc(k) + T*fxc(k); %x(n+1) = x(n) + T*x'(n)
    yp(k+1) = yc(k) + T*fyc(k);
    zp(k+1) = zc(k) + T*fzc(k);
    
    % Evaluate for correction
    fxp(k+1) = alpha*(yp(k+1)+((xp(k+1)-2*(xp(k+1)).^3)/(7)));
    fyp(k+1) = xp(k+1)-yp(k+1)+zp(k+1);
    fzp(k+1) = -(100/7)*yp(k+1);
    % Correct Trapezoidal Method
    xc(k+1) = xc(k) + (T/2)*(fxp(k+1) + fxp(k)); %x(n+1) = x(n) + (T/2) * (x'(n+1) +x'(n))
    yc(k+1) = yc(k) + (T/2)*(fyp(k+1) + fyp(k));
    zc(k+1) = zc(k) + (T/2)*(fzp(k+1) + fzp(k));
end

figure;
subplot(2,2,1)
plot3(xc,yc,zc)
title('Chuas System: 3D Plot')


subplot(2,2,2)
plot(xc, t)
ylim([0,10])
title('X vs T')

subplot(2,2,3)
plot(yc, t)
ylim([0,10])
title('Y vs T')

subplot(2,2,4)
plot(zc, t)
ylim([0,10])
title('Z vs T')