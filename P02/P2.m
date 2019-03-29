%Tyler Matthews 1/244/19
%System Simulation Problem 02

clc; %CLEAR CONSOLE 

%INITIAL VARIABLES
z = 0.5;
N = 500;
t = linspace(0,z,N);
T = z/N;

%X
A = ((13333.3/2)*T^3 + 33.3333*T^2 - T);
B = ((13333.3/2)*T^3 - 33.3333*T^2 + 2*T);
C = (-1*T);

%Y
D = ((693333/2)*T^3 + 12133.3*T^2 + 185.333*T - 3);
E = ((693333/2)*T^3 - 12133.3*T^2 - (185.333*2)*T + 3);
F = (185.333*T - 1);

%STEP RESPONSE ------------------------------------------------------------
X=ones(1,N);
Y=zeros(1,N);

Y(1) = 50*(0*A + 0*B + 0*C) - (0*F + 0*E + 0*D);
Y(2) = 50*(X(1)*A + 0*B + 0*C) - (0*F + 0*E + Y(1)*D);
Y(3) = 50*(X(2)*A + X(1)*B + 0*C) - (0*F + Y(1)*E + Y(2)*D);
% Y[4] = 50*(1*A + 2*B + 1*C) - (Y[0]*F + Y[1]*E + Y[2]*D);

%Looping Through Difference Equation
for k=4:N-1
    Y(k) = 50*(X(k-1)*A + X(k-2)*B + X(k-3)*C) - (Y(k-3)*F + Y(k-2)*E + Y(k-1)*D);
end

%Plotting
figure;
plot(t,Y);
xlabel('t');
ylabel('Y');
title('Step Response');
%END STEP RESPONSE --------------------------------------------------------


%IMPULSE RESPONSE ---------------------------------------------------------
Y=zeros(1,N);
X = zeros(1,N);
X(1) = 1/T;

Y(1) = 50*(0*A + 0*B + 0*C) - (0*F + 0*E + 0*D);
Y(2) = 50*(X(1)*A + 0*B + 0*C) - (0*F + 0*E + Y(1)*D);
Y(3) = 50*(X(2)*A + X(1)*B + 0*C) - (0*F + Y(1)*E + Y(2)*D);

%Looping Through Difference Equation
for k=4:N-1
    Y(k) = 50*(X(k-1)*A + X(k-2)*B + X(k-3)*C) - (Y(k-3)*F + Y(k-2)*E + Y(k-1)*D);
end

%Plotting
figure;
plot(t,Y);
xlabel('t');
ylabel('Y');
title('Impulse Response');
%END IMPULSE RESPONSE -----------------------------------------------------
