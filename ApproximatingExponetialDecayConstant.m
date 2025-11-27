
%The aim of this program is to find an approximation of a one-dimensional
%exponential decay model where the decay constant, k, is unknown. In order
%to mimic real world conditions, we will determine this approximation from
%data that contains noise. The generation of this data begins on line 25.
%Note that 'Kt' is the true decay constant.  

Kt = 1.5; 

%We define our time, time step, and number of steps for the problem at
%hand. The user has freedom to choose the interval of time and the time
%step, the rest will follow necessarily. 

%Let 'a' be the lower bound of the time interval and 'b', the upper bound. We
%let 'n' denote the number of subintervals induced by a step size 'h'.

a = 0; 
b = 2; 
h = 0.2; 

t = a:h:b; 
n = (b-a)/h;

%We now generate the noisy data values for 'y'. We call this set, 'Ydata'.  

Ydata = zeros(11,1);

for i = 1:n+1
    Ydata(i) = exp(-Kt*t(i))  + normrnd(0,1)/4;
    %0.1; (this was the noise used to experimentally determine a good
    %approximaiton for tau where tau is the optimal step size in the 
    %steepest descent method).
end

%In order to approximate 'k', we find the value of 'k' that will minimize the
%sum of least square deviations between 'Ydata' and the exact 'Y' which we
%lablel as 'Ysoln'. This least squares funciton will be labeled by 'J' then
%followed by a specific context, e.g. consider line 78 or line 98. 

%For the aforementioned minimization, we apply the steepest descent method 
%with tau = 0.1, a selection that was made after comparing the convergence 
%rate of other values of tau. Note that 'Yinit' is the initial condition of 
%the given differential equation. Note that k=0 is an initial guess for 'k'.

tau = 0.1;
k = 0;
Ysoln = zeros(n+1,1);
Yinit = [0 1];
maxIterations = 500;
tol = 10^(-3);

%Since I needed to store the values of 'J(kmin)' for each iteration in 
%order to observe the convergence of 'J(k)' per iteration step, I decided 
%to create a storage vector. This vector is listed as 'JkminStorage'. 
%However, the size of this vector does not correspond to the number of 
%actual approximations obtained. Hence, I introduce a counter by the name 
%of 'Iterations' in order to count the total number of steps that are 
%actually needed to obtain the acceptable approximation. The value of 
%'Iterations' is then the size of a new vector known below as 'Jplot'. 
%I then copy the nonzero elements from 'JkminStorage' to 'Jplot' in order 
%to have a plot that compares 'J(k)' to each iteration of the loop below. 

Iterations = 0;
JkminStorage = zeros(maxIterations,1);

for i = 1:maxIterations 

    Iterations = Iterations+1;

    [t,y] = ode45(@(t,y,kvalue) ODE(y,k), t, Yinit);

    Ysoln(:,1) = y(:,2);

    %Before approximating the gradient of 'J' with respect to 'k', we 
    %determine the value of 'J' given a value for 'k' for each time step. 
    %This is a crucial value since we will use the forward and backward 
    %difference formulas to obtain the approximation of the derivative of 
    %'J' at 'k'.

    Jk = (1/2)*(Ysoln(1,1)-Ydata(1))^2;

    for p = 2:n+1
        Jk = Jk + (1/2)*(Ysoln(p,1)-Ydata(p))^2; 
    end

    %To compute the gradient of 'J' at 'k' using a finite difference
    %method, we need at least one other value for 'k'. Hence, we use some
    %step size called 'dk' to get a nearby 'k' that we call 'kPlus' from 
    %which we can approximate the derivative of 'J' at 'k'. 

    dk =0.1;
    kPlus = k +dk; 

    [t,y] = ode45(@(t,y,kvalue) ODE(y,kPlus), t, Yinit);
    Ysoln(:,1) = y(:,2);


    %We now find the corresponding value of 'J' which we call 'JkPlus'.

    JkPlus = (1/2)*(Ysoln(1,1)-Ydata(1))^2;

    for p = 2:n+1
        JkPlus = JkPlus + (1/2)*(Ysoln(p,1)-Ydata(p))^2; 
    end

    %Using the forward difference formula, we obtain the derivative of 'J' 
    %at some 'k', we call this derivative 'JkPrime'. Note that this 
    %derivative is the gradient of 'J' at 'k'. 

    JkPrime = (JkPlus-Jk)/dk;

    %The next line finds the next approximation for the 'k' that minimizes
    %the objection function. 
   
    kmin = k - tau*JkPrime;

    %We now use an ODE solver to find the value of 'J' at 'kmin'.

    [t,y] = ode45(@(t,y,kvalue) ODE(y,kmin), t, Yinit);
    Ysoln(:,1) = y(:,2);

    Jkmin = (1/2)*(Ysoln(1,1)-Ydata(1))^2;

    for p = 2:n+1
        Jkmin = Jkmin + (1/2)*(Ysoln(p,1)-Ydata(p))^2; 
    end

    JkminStorage(i,1) = Jkmin;

    %We check if the previous approximation and the current approximation
    %are sufficiently close to be considered just about the same point. 

    if abs(kmin - k) < tol
        Iterations
        kmin
        break
    end

    k = kmin;
end

%We now compare our approximations to the actual solution by way of graphs. 
%The titles below describe these graphs. We make a new vector titled
%'iterations' that will become the x-axis on the graph that compares the
%value of 'J(k)' for each iteration towards convergence 'J(kmin)'. 

   Jplot = zeros(Iterations, 1);
   iterations = ones(Iterations, 1);

   for i = 1:Iterations
       Jplot(i) = JkminStorage(i);
   end

   for i = 2:Iterations
       iterations(i) = iterations(i-1)+1;
   end

   %In the following, we determine the fitted model that is generated by
   %'kmin' that was found using the steepest descent method. 

   [t,y] = ode45(@(t,y,kvalue) ODE(y,kmin), t, Yinit);
    Ykmin(:,1) = y(:,2);


    figure(1)
    plot(t, Ydata)
    hold on 
    plot(t, Ykmin)
    legend('Data', 'Fitted Model Using $k_{approx}$', 'interpreter', 'latex');
    xlabel('$t$', 'interpreter', 'latex')
    ylabel('$y(t;k)$', 'interpreter', 'latex')
    title('Comparing The Data and The Fitted Model')
    hold off


    figure(2)
    plot(iterations, Jplot)
    hold on
    xlabel('Iteration Number')
    ylabel('$J(k)$', 'interpreter', 'latex')
    title('Convergence of $J(k)$ to Minimum as Iterations Increase', 'interpreter', 'latex')
    hold off
    
function dydt = ODE(y,k)

    dydt = -k*y;

end