%Tyler Matthews
%P15

clc; close all; clear all; %Clear console, clear workspace, and close figures

H = 100;

%% Part A
R0 = 0.5;

% Case 1
a1 = 0.3;
b1 = 0.2;
c1 = 0.5;
d1 = 0.03302;
e1 = 0.01;
m01 = R0/((a1^2)*b1*c1*(1/e1)*(1/d1))

% Case 2
a2 = 0.1;
b2 = 0.03;
c2 = 0.275;
d2 = 0.03304;
e2 = 0.0035;
m02 = R0/((a2^2)*b2*c2*(1/e2)*(1/d2))

% Case 3
a3 = 0.5;
b3 = 0.4;
c3 = 0.4;
d3 = 0.1;
e3 = 0.05;
m03 = R0/((a3^2)*b3*c3*(1/e3)*(1/d3))

%% Part B
T = 0.0001;
t0plot = 10;
tfinal = 250;
t = 0:T:tfinal;
k0 = 1+floor(t0plot/T);
N = length(t);

% For m < m0
% Case 1
m1 = 0.001;
V1 = m1*H;

y1 = ones(1,N);
i1 = ones(1,N);
y11 = y1;
i21 = i1;
y12 = y1;
i22 = i1;
f11 = ones(1,N);
f21 = ones(1,N);
f12 = f11;
f22 = f21;

for k=1:N-1
    f11(k+1) = a1*b1*i1(k)*((H-y1(k))/H) - e1*y1(k);
    f21(k+1) = a1*c1*(V1-i1(k))*(y1(k)/H) - d1*i1(k);
    
    y11(k+1) = y1(k) + 0.5*T*f11(k+1);
    i21(k+1) = i1(k) + 0.5*T*f21(k+1);
    
    f12(k+1) = a1*b1*i21(k+1)*((H-y11(k+1))/H) - e1*y11(k+1);
    f22(k+1) = a1*c1*(V1-i21(k+1))*(y11(k+1)/H) - d1*i21(k+1);
    
    y1(k+1) = y1(k) + T*f12(k+1);
    i1(k+1) = i1(k) + T*f22(k+1);
end

figure;
subplot(2,3,1)
plot(t,y1)
xlabel('t')
ylabel('Y')
title('Infected Humans m=0.001')

subplot(2,3,4)
plot(t,i1)
xlabel('t')
ylabel('I')
title('Infected Mosquitos m=0.001')

% Case 2
m2 = 0.007;
V2 = m2*H;

y2 = ones(1,N);
i2 = ones(1,N);
y31 = y2;
i41 = i2;
f31 = ones(1,N);
f41 = ones(1,N);
f32 = f31;
f42 = f41;

for k=1:N-1
    f31(k+1) = a1*b1*i2(k)*((H-y2(k))/H) - e1*y2(k);
    f41(k+1) = a1*c1*(V2-i2(k))*(y2(k)/H) - d1*i2(k);
    
    y31(k+1) = y2(k) + 0.5*T*f31(k+1);
    i41(k+1) = i2(k) + 0.5*T*f41(k+1);
    
    f32(k+1) = a1*b1*i41(k+1)*((H-y31(k+1))/H) - e1*y31(k+1);
    f42(k+1) = a1*c1*(V2-i41(k+1))*(y31(k+1)/H) - d1*i41(k+1);
    
    y2(k+1) = y2(k) + T*f32(k+1);
    i2(k+1) = i2(k) + T*f42(k+1);
end

subplot(2,3,2)
plot(t,y2)
xlabel('t')
ylabel('Y')
title('Infected Humans m=0.007')

subplot(2,3,5)
plot(t,i2)
xlabel('t')
ylabel('I')
title('Infected Mosquitos m=0.007')

% Case 3
m3 = 0.006;
V3 = m3*H;

y3 = ones(1,N);
i3 = ones(1,N);
y51 = y3;
i61 = i3;
f51 = ones(1,N);
f61 = ones(1,N);
f52 = f51;
f62 = f61;

for k=1:N-1
    f51(k+1) = a1*b1*i3(k)*((H-y3(k))/H) - e1*y3(k);
    f61(k+1) = a1*c1*(V3-i3(k))*(y3(k)/H) - d1*i3(k);
    
    y51(k+1) = y3(k) + 0.5*T*f51(k+1);
    i61(k+1) = i3(k) + 0.5*T*f61(k+1);
    
    f52(k+1) = a1*b1*i61(k+1)*((H-y51(k+1))/H) - e1*y51(k+1);
    f62(k+1) = a1*c1*(V3-i61(k+1))*(y51(k+1)/H) - d1*i61(k+1);
    
    y3(k+1) = y3(k) + T*f52(k+1);
    i3(k+1) = i3(k) + T*f62(k+1);
end

subplot(2,3,3)
plot(t,y3)
xlabel('t')
ylabel('Y')
title('Infected Humans m=0.006')

subplot(2,3,6)
plot(t,i3)
xlabel('t')
ylabel('I')
title('Infected Mosquitos m=0.006')

% For m > m0
% Case 1
m11 = 0.02;
V11 = m11*H;

y1 = ones(1,N);
i1 = ones(1,N);
y11 = y1;
i21 = i1;
y12 = y1;
i22 = i1;
f11 = ones(1,N);
f21 = ones(1,N);
f12 = f11;
f22 = f21;

for k=1:N-1
    f11(k+1) = a1*b1*i1(k)*((H-y1(k))/H) - e1*y1(k);
    f21(k+1) = a1*c1*(V11-i1(k))*(y1(k)/H) - d1*i1(k);
    
    y11(k+1) = y1(k) + 0.5*T*f11(k+1);
    i21(k+1) = i1(k) + 0.5*T*f21(k+1);
    
    f12(k+1) = a1*b1*i21(k+1)*((H-y11(k+1))/H) - e1*y11(k+1);
    f22(k+1) = a1*c1*(V11-i21(k+1))*(y11(k+1)/H) - d1*i21(k+1);
    
    y1(k+1) = y1(k) + T*f12(k+1);
    i1(k+1) = i1(k) + T*f22(k+1);
end

figure;
subplot(2,3,1)
plot(t,y1)
xlabel('t')
ylabel('Y')
title('Infected Humans m=0.02')

subplot(2,3,4)
plot(t,i1)
xlabel('t')
ylabel('I')
title('Infected Mosquitos m=0.02')

% Case 2
m22 = 0.8;
V22 = m22*H;

y2 = ones(1,N);
i2 = ones(1,N);
y31 = y2;
i41 = i2;
f31 = ones(1,N);
f41 = ones(1,N);
f32 = f31;
f42 = f41;

for k=1:N-1
    f31(k+1) = a1*b1*i2(k)*((H-y2(k))/H) - e1*y2(k);
    f41(k+1) = a1*c1*(V22-i2(k))*(y2(k)/H) - d1*i2(k);
    
    y31(k+1) = y2(k) + 0.5*T*f31(k+1);
    i41(k+1) = i2(k) + 0.5*T*f41(k+1);
    
    f32(k+1) = a1*b1*i41(k+1)*((H-y31(k+1))/H) - e1*y31(k+1);
    f42(k+1) = a1*c1*(V22-i41(k+1))*(y31(k+1)/H) - d1*i41(k+1);
    
    y2(k+1) = y2(k) + T*f32(k+1);
    i2(k+1) = i2(k) + T*f42(k+1);
end

subplot(2,3,2)
plot(t,y2)
xlabel('t')
ylabel('Y')
title('Infected Humans m=0.8')

subplot(2,3,5)
plot(t,i2)
xlabel('t')
ylabel('I')
title('Infected Mosquitos m=0.8')

% Case 3
m33 = 0.07;
V33 = m33*H;

y3 = ones(1,N);
i3 = ones(1,N);
y51 = y3;
i61 = i3;
f51 = ones(1,N);
f61 = ones(1,N);
f52 = f51;
f62 = f61;

for k=1:N-1
    f51(k+1) = a1*b1*i3(k)*((H-y3(k))/H) - e1*y3(k);
    f61(k+1) = a1*c1*(V33-i3(k))*(y3(k)/H) - d1*i3(k);
    
    y51(k+1) = y3(k) + 0.5*T*f51(k+1);
    i61(k+1) = i3(k) + 0.5*T*f61(k+1);
    
    f52(k+1) = a1*b1*i61(k+1)*((H-y51(k+1))/H) - e1*y51(k+1);
    f62(k+1) = a1*c1*(V33-i61(k+1))*(y51(k+1)/H) - d1*i61(k+1);
    
    y3(k+1) = y3(k) + T*f52(k+1);
    i3(k+1) = i3(k) + T*f62(k+1);
end

subplot(2,3,3)
plot(t,y3)
xlabel('t')
ylabel('Y')
title('Infected Humans m=0.07')

subplot(2,3,6)
plot(t,i3)
xlabel('t')
ylabel('I')
title('Infected Mosquitos m=0.07')

% For m >> m0
% Case 1
m11 = 0.1;
V11 = m11*H;

y1 = ones(1,N);
i1 = ones(1,N);
y11 = y1;
i21 = i1;
y12 = y1;
i22 = i1;
f11 = ones(1,N);
f21 = ones(1,N);
f12 = f11;
f22 = f21;

for k=1:N-1
    f11(k+1) = a1*b1*i1(k)*((H-y1(k))/H) - e1*y1(k);
    f21(k+1) = a1*c1*(V11-i1(k))*(y1(k)/H) - d1*i1(k);
    
    y11(k+1) = y1(k) + 0.5*T*f11(k+1);
    i21(k+1) = i1(k) + 0.5*T*f21(k+1);
    
    f12(k+1) = a1*b1*i21(k+1)*((H-y11(k+1))/H) - e1*y11(k+1);
    f22(k+1) = a1*c1*(V11-i21(k+1))*(y11(k+1)/H) - d1*i21(k+1);
    
    y1(k+1) = y1(k) + T*f12(k+1);
    i1(k+1) = i1(k) + T*f22(k+1);
end

figure;
subplot(2,3,1)
plot(t,y1)
xlabel('t')
ylabel('Y')
title('Infected Humans m=0.1')

subplot(2,3,4)
plot(t,i1)
xlabel('t')
ylabel('I')
title('Infected Mosquitos m=0.1')

% Case 2
m22 = 7;
V22 = m22*H;

y2 = ones(1,N);
i2 = ones(1,N);
y31 = y2;
i41 = i2;
f31 = ones(1,N);
f41 = ones(1,N);
f32 = f31;
f42 = f41;

for k=1:N-1
    f31(k+1) = a1*b1*i2(k)*((H-y2(k))/H) - e1*y2(k);
    f41(k+1) = a1*c1*(V22-i2(k))*(y2(k)/H) - d1*i2(k);
    
    y31(k+1) = y2(k) + 0.5*T*f31(k+1);
    i41(k+1) = i2(k) + 0.5*T*f41(k+1);
    
    f32(k+1) = a1*b1*i41(k+1)*((H-y31(k+1))/H) - e1*y31(k+1);
    f42(k+1) = a1*c1*(V22-i41(k+1))*(y31(k+1)/H) - d1*i41(k+1);
    
    y2(k+1) = y2(k) + T*f32(k+1);
    i2(k+1) = i2(k) + T*f42(k+1);
end

subplot(2,3,2)
plot(t,y2)
xlabel('t')
ylabel('Y')
title('Infected Humans m=7')

subplot(2,3,5)
plot(t,i2)
xlabel('t')
ylabel('I')
title('Infected Mosquitos m=7')

% Case 3
m33 = 0.6;
V33 = m33*H;

y3 = ones(1,N);
i3 = ones(1,N);
y51 = y3;
i61 = i3;
f51 = ones(1,N);
f61 = ones(1,N);
f52 = f51;
f62 = f61;

for k=1:N-1
    f51(k+1) = a1*b1*i3(k)*((H-y3(k))/H) - e1*y3(k);
    f61(k+1) = a1*c1*(V33-i3(k))*(y3(k)/H) - d1*i3(k);
    
    y51(k+1) = y3(k) + 0.5*T*f51(k+1);
    i61(k+1) = i3(k) + 0.5*T*f61(k+1);
    
    f52(k+1) = a1*b1*i61(k+1)*((H-y51(k+1))/H) - e1*y51(k+1);
    f62(k+1) = a1*c1*(V33-i61(k+1))*(y51(k+1)/H) - d1*i61(k+1);
    
    y3(k+1) = y3(k) + T*f52(k+1);
    i3(k+1) = i3(k) + T*f62(k+1);
end

subplot(2,3,3)
plot(t,y3)
xlabel('t')
ylabel('Y')
title('Infected Humans m=0.6')

subplot(2,3,6)
plot(t,i3)
xlabel('t')
ylabel('I')
title('Infected Mosquitos m=0.6')

%% Part C
disp('Case 1: The infection is very difficult to transmit beacuse all of the graphs in Figure 1 decay.')

disp('Case 2: The infection is relatively easy to transmit because it only decays when m<m0 as can be seen in Figure 2.')

disp('Case 3: The infection can be transmitted easily because none of the graphs in Figure 3 decay.')
