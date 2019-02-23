%System Simulation Problem 3
%Tyler Matthews 2/6/2019

%PART A
%Initial Coniditions
z = 0.05;
N = 5000;
t = linspace(0,z,N);
T = z/N;

%Coefficients
A = (4/3)*1e7;
B = (1/(T^2));
C = ((-2/(T^2)) - (250/T));
D = ((1/(T^2)) + (250/T) + (3.33e7));

%Input & Ouput Arrays
Vi = 12 * ones(1,N);
Vo = zeros(1,N);

%First two places in output array
Vo(1) = (Vi(1)*A) / (D);
Vo(2) = (Vi(2)*A - Vo(1)*C) / (D);

%Calculating the rest of the output array
for k=3:N-1
    Vo(k) = (Vi(k)*A - Vo(k-1)*C - Vo(k-2)*B) / (D);
end

%Plotting
plot(t, Vo);
title('Problem 3A');
xlabel('t');
ylabel('Vo');
figure;

%PART B
%Initial Coniditions
z = 0.05;
N = 5000;
t = linspace(0,z,N);
T = z/N;

%Coefficients
A = 4e8;
B = (1/(T^2));
C = ((-2/(T^2)) - (250/T));
D = ((1/(T^2)) + (250/T) + (3.33e7));

%Input & Ouput Arrays
Vi = 0.4 * ones(1,N);
Vo = zeros(1,N);

%First two places in output array
Vo(1) = (Vi(1)*A) / (D);
Vo(2) = (Vi(2)*A - Vo(1)*C) / (D);

%Calculating the rest of the output array
for k=3:N-1
    Vo(k) = (Vi(k)*A - Vo(k-1)*C - Vo(k-2)*B) / (D);
end

%Plotting
plot(t, Vo);
title('Problem 3B');
xlabel('t');
ylabel('Vo');