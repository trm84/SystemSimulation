%Tyler Matthews P10

%% Part A
clc; close all; %Clear Console, close figures

num = [0 0 0 0 0.0850];
den = [1 0.4174 1.0871 0.2805 0.1512];

poles = roots(den)

%% Part B
figure;
G = tf(num, den)
bode(G)

%% PART C
abDen = [1 -1 0];
abNum = [0 3 -1];

Phi = tf(abDen, abNum) %(sigma / roe) : (row - l*sigma)
newNum = [3 -2 1];
badPoints = roots(newNum)
magnitude = abs(badPoints)

figure;
Nt=21;
Nr=12;

theta=linspace(0,2*pi,1001); 
rho=linspace(0.5774,1,1001); 
tvec=linspace(0,2*pi,Nt);
rvec=linspace(0.5774,1,Nr);

for k=1:length(rvec)
 z=rvec(k)*exp(i*theta);
 w=(z.^2-z)./(3.*z - 1);
 
 hold on
 plot(real(w), imag(w))
 hold off
end

for k=1:length(tvec)-1
 z=rho*exp(i*tvec(k));
 w=(z.^2-z)./(3.*z - 1);
 hold on
 plot(real(w), imag(w))
 hold off
end
title('AB2 Stability Region')

figure;
Nt=21;
Nr=12;

theta=linspace(0,2*pi,1001); 
rho=linspace(0, 0.5774,1001); 
tvec=linspace(0,2*pi,Nt);
rvec=linspace(0, 0.5774,Nr);

for k=1:length(rvec)
 z=rvec(k)*exp(i*theta);
 w=(z.^2-z)./(3.*z - 1);
 
 hold on
 plot(real(w), imag(w))
 hold off
end

for k=1:length(tvec)-1
 z=rho*exp(i*tvec(k));
 w=(z.^2-z)./(3.*z - 1);
 hold on
 plot(real(w), imag(w))
 hold off
end

axis([-1 1 -1 1])

%% PART D
figure;

lamda = eig(G)
T = linspace(0,1,1001);

hold on
plot(real(lamda(1)*T), imag(lamda(1)*T))
plot(real(lamda(2)*T), imag(lamda(2)*T))
plot(real(lamda(3)*T), imag(lamda(3)*T))
plot(real(lamda(4)*T), imag(lamda(4)*T))
hold off

Nt=21;
Nr=12;

theta=linspace(0,2*pi,1001); 
rho=linspace(0.5774,1,1001); 
tvec=linspace(0,2*pi,Nt);
rvec=linspace(0.5774,1,Nr);

for k=1:length(rvec)
 z=rvec(k)*exp(i*theta);
 w=(z.^2-z)./(3.*z - 1);
 
 hold on
 plot(real(w), imag(w))
 hold off
end

for k=1:length(tvec)-1
 z=rho*exp(i*tvec(k));
 w=(z.^2-z)./(3.*z - 1);
 hold on
 plot(real(w), imag(w))
 hold off
end
axis([-0.5 0.1 -0.6 0.6])
T_stable=0.08 %Stable
T_rel_unstable=0.6 %Relatively Unstable
T_unstable= 1 %Unstable

hold on
for k = 1 : 4
    plot(T_stable*lamda(k), 'x')
    plot(T_rel_unstable*lamda(k), 'o')
    plot(T_unstable*lamda(k), '*') %Off page
end
hold off

%% PART E
poles_stable = roots(abDen - abNum*T_stable)
poles_rel_unstable = roots(abDen - abNum*T_rel_unstable)
poles_unstable = roots(abDen - abNum*T_unstable)

%% PART F

[A, B, C, D] = tf2ss(num, den)

N = 10000;
t = linspace(0,10,N);
u = ones(1,N);
fx1 = zeros(1,N);
fx2 = zeros(1,N);
fx3 = zeros(1,N);
fx4 = zeros(1,N);

x1 = zeros(1,N);
x2 = zeros(1,N);
x3 = zeros(1,N);
x4 = zeros(1,N);
y = zeros(1,N);

T=T_stable;
for k = 1:N-1
    fx1(k+1) = -0.4174*x1(k+1)-1.0871*x2(k+1)-0.2805*x3(k+1)-0.1512*x4(k+1)+u(k+1);
    fx2(k+1) = x1(k+1);
    fx3(k+1) = x2(k+1);
    fx4(k+1) = x3(k+1);

    x1(k+2) = x1(k+1) + (T/2) * (3*fx1(k+1) - fx1(k));
    x2(k+2) = x2(k+1) + (T/2) * (3*fx2(k+1) - fx2(k));
    x3(k+2) = x3(k+1) + (T/2) * (3*fx3(k+1) - fx3(k));
    x4(k+2) = x4(k+1) + (T/2) * (3*fx4(k+1) - fx4(k));
    
    y(k) = 0.085*x3(k);
end

figure
plot(t,y)
xlim([0 2])
title('Relatively Stable T')




N = 10000;
t = linspace(0,10,N);
u = ones(1,N);
fx1 = zeros(1,N);
fx2 = zeros(1,N);
fx3 = zeros(1,N);
fx4 = zeros(1,N);

x1 = zeros(1,N);
x2 = zeros(1,N);
x3 = zeros(1,N);
x4 = zeros(1,N);
y = zeros(1,N);

T=T_rel_unstable;
for k = 1:N-1
    fx1(k+1) = -0.4174*x1(k+1)-1.0871*x2(k+1)-0.2805*x3(k+1)-0.1512*x4(k+1)+u(k+1);
    fx2(k+1) = x1(k+1);
    fx3(k+1) = x2(k+1);
    fx4(k+1) = x3(k+1);

    x1(k+2) = x1(k+1) + (T/2) * (3*fx1(k+1) - fx1(k));
    x2(k+2) = x2(k+1) + (T/2) * (3*fx2(k+1) - fx2(k));
    x3(k+2) = x3(k+1) + (T/2) * (3*fx3(k+1) - fx3(k));
    x4(k+2) = x4(k+1) + (T/2) * (3*fx4(k+1) - fx4(k));
    
    y(k) = 0.085*x3(k);
end

figure
plot(t,y)
xlim([0 2])
title('Relatively Unstable T')

N = 10000;
t = linspace(0,10,N);
u = ones(1,N);
fx1 = zeros(1,N);
fx2 = zeros(1,N);
fx3 = zeros(1,N);
fx4 = zeros(1,N);

x1 = zeros(1,N);
x2 = zeros(1,N);
x3 = zeros(1,N);
x4 = zeros(1,N);
y = zeros(1,N);

T=T_unstable;
for k = 1:N-1
    fx1(k+1) = -0.4174*x1(k+1)-1.0871*x2(k+1)-0.2805*x3(k+1)-0.1512*x4(k+1)+u(k+1);
    fx2(k+1) = x1(k+1);
    fx3(k+1) = x2(k+1);
    fx4(k+1) = x3(k+1);

    x1(k+2) = x1(k+1) + (T/2) * (3*fx1(k+1) - fx1(k));
    x2(k+2) = x2(k+1) + (T/2) * (3*fx2(k+1) - fx2(k));
    x3(k+2) = x3(k+1) + (T/2) * (3*fx3(k+1) - fx3(k));
    x4(k+2) = x4(k+1) + (T/2) * (3*fx4(k+1) - fx4(k));
    
    y(k) = 0.085*x3(k);
end

figure
plot(t,y)
xlim([0 2])
title('Unstable T --> NOT REQUIRED')