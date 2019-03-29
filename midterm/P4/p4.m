%Tyler Matthews
%System Simluation Midterm P4
clc; close all; %Clear console and close figures


%% PART A
disp('Corrector is fourth order accurate because up to C3 in Lambert''s equations are equal to 0');

%% PART C
Nt=21;
Nr=12;

num = [0.46 0.29 -0.32]
den = [1 -1.56 0.56]

Hp = tf(num, den)
zeros = roots(num)
poles = roots(den)

Phi = tf(den, num) %(sigma / roe) : (row - l*sigma)

newNum = [10076 -11552 3368]
badPoints = roots(newNum)
magnitude = abs(badPoints)

theta=linspace(0,2*pi,1001); 
rho=linspace(0,0.5686,1001); 
tvec=linspace(0,2*pi,Nt);
rvec=linspace(0,0.5686,Nr);

temp = (roots(den - num*0.5764));
mag = abs(temp)
ang = angle(temp)

for k=1:length(rvec)
 z=rvec(k)*exp(i*theta);
 w = (z.^2-z.*1.56 + 0.56)./(z.^2.*0.46+z.*0.29-0.32);
 hold on
 plot(real(w), imag(w))
 hold off
end

for k=1:length(tvec)-1
 z=rho*exp(i*tvec(k));
 w = (z.^2-z.*1.56 + 0.56)./(z.^2.*0.46+z.*0.29-0.32);
 hold on
 plot(real(w), imag(w))
 hold off
end

grid on
axis([-6 2.5 -2.5 2.5])
title('Primary Domain')

%TESTING FOR INTERSECTION POINT
% for N=1:10
%     temp =  0.568 + N*0.0001
%     val = sprintf('N = %0.5f',temp);
%     z = (temp) * exp(i*theta);
%     w = (z.^2-z.*1.56 + 0.56)./(z.^2.*0.46+z.*0.29-0.32);
%     plot(real(w), imag(w));
%     title(val);
%     disp(val);
%     disp(w(1));
%     disp(w(2));
%     pause;
% end

%% PART D
%Getting an actual value for T using the state-space representation in
%problem 03

A = [-4.7, -1.55, -0.55; 0.3, -2.75, -0.35; 1.1, 1.85, -2.55]
B = [1; 0; -1]
C = [2, 1, 1]
D = [0]

lamda = eig(A)

figure;

for k=1:length(rvec)
 z=rvec(k)*exp(i*theta);
 w = (z.^2-z.*1.56 + 0.56)./(z.^2.*0.46+z.*0.29-0.32);
 hold on
 plot(real(w), imag(w))
 hold off
end

for k=1:length(tvec)-1
 z=rho*exp(i*tvec(k));
 w = (z.^2-z.*1.56 + 0.56)./(z.^2.*0.46+z.*0.29-0.32);
 hold on
 plot(real(w), imag(w))
 hold off
end

axis([-6 2.5 -2.5 2.5])
title('Primary Domain with Lamda-T Products')
T = linspace(0,2,1001);

hold on
plot(real(lamda(1)*T), imag(lamda(1)*T))
plot(real(lamda(2)*T), imag(lamda(2)*T))
plot(real(lamda(3)*T), imag(lamda(3)*T))
hold off

disp('Stable, and accurate, values for T would be located inside of the green incompassing line and to the right of the plot''s origin (-1.75, 0)')
disp('The best points would be between -1 & -1.75 on the x-axis and between 1 & -1 on the y-axis')
disp('Using lamda-T products from Problem 3, this would give us T = 0.23')

%% PART E
figure;

for k=1:length(rvec)
 z=rvec(k)*exp(i*theta);
 w = (z.^2-z.*1.56 + 0.56)./(z.^2.*0.46+z.*0.29-0.32);
 hold on
 plot(real(w), imag(w), 'r')
 hold off
end

for k=1:length(tvec)-1
 z=rho*exp(i*tvec(k));
 w = (z.^2-z.*1.56 + 0.56)./(z.^2.*0.46+z.*0.29-0.32);
 hold on
 plot(real(w), imag(w), 'r')
 hold off
end

theta=linspace(0,2*pi,1001); 
rho=linspace(0,0.5256,1001); 
tvec=linspace(0,2*pi,Nt);
rvec=linspace(0,0.5256,Nr);


for k=1:length(rvec)
 z=rvec(k)*exp(i*theta);
 w=(z.^2-z.*1.45 + 0.45)./(z.*1.27-0.73);
 
 hold on
 plot(real(w), imag(w), 'b')
 hold off
end

for k=1:length(tvec)-1
 z=rho*exp(i*tvec(k));
 w=(z.^2-z.*1.45 + 0.45)./(z.*1.27-0.73);
 hold on
 plot(real(w), imag(w), 'b')
 hold off
end


grid on
axis([-6 2.5 -2.5 2.5])
title('Corrector and Predictor Stability Plots')

disp(' ')
disp('Stability plot for the predictor is in blue and the stability plot for corrector is in red')
disp('You can see that the plots are similiar, but the corrector plot is much larger and shifted to the left')
