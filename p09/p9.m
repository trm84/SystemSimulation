%Tyler Matthews
%System Simulation Problem 9

clc; %clear console and close figures

num = [5.269 3.076 -0.731];
den = [1 6.81 -9 1.193];

integrator_zeros = roots(num)
integrator_poles = roots(den)

%derviative of den/num = newNum/newDen
newNum = [5.267 6.15 66.138 -22.517 2.91]
newDen = [27.742 32.39 1.755 -4.49 0.534]

badPoints = roots(newNum)
magnitude = abs(badPoints)

Nt=21;
Nr=12;

theta=linspace(0,2*pi,1001); 
rho=linspace(0.2060,3.6076,1001); 
tvec=linspace(0,2*pi,Nt);
rvec=linspace(0.2060,3.6076,Nr);

figure(2);


for k=1:length(rvec)
 z=rvec(k)*exp(i*theta);
 w=(z.^3 + 6.81*z.^2 - 9*z + 1.193)./(5.269*z.^2 + 3.076*z - 0.731);
 
 hold on
 plot(real(w), imag(w))
 hold off
end

for k=1:length(tvec)-1
 z=rho*exp(i*tvec(k));
 w=(z.^3 + 6.81*z.^2 - 9*z + 1.193)./(5.269*z.^2 + 3.076*z - 0.731);
 
 hold on
 plot(real(w), imag(w))
 hold off
end
axis([-5 5 -10 10])