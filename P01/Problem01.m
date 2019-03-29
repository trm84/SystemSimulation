%Tyler Matthews 1/16/19
%System Simulation Problem 01

N = 201;
t = linspace(0,10,N);

%Given Values = 0.8, 1.35, 2.75, 3.2, 3.52, 4,
%Chosen Values = 2.9, 3, 3.7 4.1
a = [0.8, 1.35, 2.75, 3.2, 3.52, 4, 2.9, 3, 3.7, 4.1];

%Loops through all values of a
for i=1:10
    %Initial Conditions
    x = zeros(1,N);
    x(1) = 0.11;

    %Difference Equation
    for k=1:N-1
       x(k+1) = a(i) * (1-x(k))*x(k);
    end
    
    %Plotting
    figure;
    plot (t,x);
    xlabel('t');
    ylabel('x');
    plot_title = sprintf('%s %0.2f','a = ',a(i));
    title(plot_title);
end

