%Tyler Matthews
%2/23/2019

R = 400;
L = 3e-3;
C = 10e-6;
S = 0.4;
Vin = 12;


A = [0, -1/L;
     1/C, -1/(R*C)];
C = [0, 1];
D = [0];

%PART A
%Vout / Vin
B = [S/L; 0]
[NUM,DEN] = ss2tf(A,B,C,D);
TF = tf(NUM,DEN)

%Vout / S
B = [Vin/L; 0]
[NUM,DEN] = ss2tf(A,B,C,D);
TF = tf(NUM,DEN)

POLES = pole(TF) %Only run this once because the poles do not change with a change in B
EIGENVALUES = eig(A)

%PART B
B = [1/L; 0];

T = 1/55000;
z = 1000;
t = 0 : T/z : 0.015;
N = length(t);

x = zeros(2, N);
y = zeros(1, N);
f = zeros(2, N);

%Setting Duty Cycle
%LOW
j = mod(0:N-1, 100);
temp = find(j>40);
j(temp) = 0;
%HIGH
temp = find(j);
j(temp) = 1;

%Input
Vin = 12 * ones(1,N);
u = Vin.*j;

x(1) = 0; %Inital Condition

for k=1:N-1
    f(:,k) = A*x(:,k) + B*u(k);
    x(:,k+1) = x(:,k) + (T/1000)*f(:,k);
    y(k+1) = C*x(:,k+1);
end

plot(t,y);
title('Buck Converter');
xlabel('t');
ylabel('v');
