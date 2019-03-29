%Tyler Matthews
%System Simluation Midterm P3
close; clc; close all;

%% E = 6
T = 0.0001;
t = 0:T:250;
N = length(t);

fx1 = zeros(1,N);
fx2 = zeros(1,N);
fx3 = zeros(1,N);
x1 = zeros(1,N);
x2 = zeros(1,N);
x3 = zeros(1,N);

e = 6;
w = 2*pi/10;

for k = 1:N-1
    u = 1.2*cos(w*t(k));

    fx1(k+1) = x2(k+1);
    fx2(k+1) = e*(1-x1(k+1)^2)*x2(k+1)-x1(k+1) + u;
    fx3(k+1) = w^2*u;
    
    x1(k+2) = x1(k+1) + (T/2)*(3*fx1(k+1) - fx1(k));
    x2(k+2) = x2(k+1) + (T/2)*(3*fx2(k+1) - fx2(k));
    x3(k+2) = x3(k+1) + (T/2)*(3*fx3(k+1) - fx3(k));
end
figure
plot(x1(10000:N),x2(10000:N)) %Remove Transient
title('e = 6')
disp('')


%% E = 8.53
T = 0.0001;
t = 0:T:250;
N = length(t);

fx1 = zeros(1,N);
fx2 = zeros(1,N);
fx3 = zeros(1,N);
x1 = zeros(1,N);
x2 = zeros(1,N);
x3 = zeros(1,N);

e = 8.53;
w = 2*pi/10;

for k = 1:N-1
    u = 1.2*cos(w*t(k));

    fx1(k+1) = x2(k+1);
    fx2(k+1) = e*(1-x1(k+1)^2)*x2(k+1)-x1(k+1) + u;
    fx3(k+1) = w^2*u;
    
    x1(k+2) = x1(k+1) + (T/2)*(3*fx1(k+1) - fx1(k));
    x2(k+2) = x2(k+1) + (T/2)*(3*fx2(k+1) - fx2(k));
    x3(k+2) = x3(k+1) + (T/2)*(3*fx3(k+1) - fx3(k));
end
figure
plot(x1(10000:N),x2(10000:N)) %Remove Transient
title('e = 8.53')
disp('')

%% E = 10
T = 0.0001;
t = 0:T:250;
N = length(t);

fx1 = zeros(1,N);
fx2 = zeros(1,N);
fx3 = zeros(1,N);
x1 = zeros(1,N);
x2 = zeros(1,N);
x3 = zeros(1,N);

e = 10;
w = 2*pi/10;

for k = 1:N-1
    u = 1.2*cos(w*t(k));

    fx1(k+1) = x2(k+1);
    fx2(k+1) = e*(1-x1(k+1)^2)*x2(k+1)-x1(k+1) + u;
    fx3(k+1) = w^2*u;
    
    x1(k+2) = x1(k+1) + (T/2)*(3*fx1(k+1) - fx1(k));
    x2(k+2) = x2(k+1) + (T/2)*(3*fx2(k+1) - fx2(k));
    x3(k+2) = x3(k+1) + (T/2)*(3*fx3(k+1) - fx3(k));
end
figure
plot(x1(10000:N),x2(10000:N)) %Remove Transient
title('e = 10')

%% Characterize
disp('All three plots are oscillatory in nature')
disp('Once reaching the limit cycle they continue on in the same pattern forever')
disp('E=6 is the slowest moving plot, where it only wraps around a few times')
disp('E = 8.53 is the fastest plot, where it wraps around many times. Also, there are multiple spots where it loops around itself')
disp('E=10 is faster than E=6 and slower than E=8.53. It also does not wrap around itself like it does in E=8.53')