function  FDMWaveDirichelet

%This code finds an approximation to the solution of a one-dimensional wave 
%equation. In LaTeX notation, this equation is $u_{tt} = 4u_{xx}$. The
%boundary conditions can be inferred by the lines 42-53. After obtaining
%the approximate solution, we plot a graph of this solution with respect to
%both time and space. 

h = 0.05;
k = 0.025;

%Note that the 'k' we choose assures convergence of the algorithm below to
%an approximate solution. This 'k' was inspired by the term that we label
%as 'labmda' which comes from the finit difference methods used to 
%approximate the derivatives in this problem. 

lambda = (2*k)/h;

n = 1/h +1;
m = 1/k + 1; 

n = int64(n);
m = int64(m);

x = 0:h:1;
t = 0:k:1; 


%Since isolating the unknowns in each approximation of the derivaties was
%fairly straight forward, we do not make use of a coefficient matrix to
%find an approximation of the solution.

Uapprox = zeros(n,m);
Uexact = zeros(n,m);

for i = 1:n
    for j = 1:m
        Uexact(i,j) = (x(i)+2*t(j))^5;
    end
end



for i = 1:n
    Uapprox(i,1) = (x(i))^5;
end

for j = 1:m
    Uapprox(1,j) = 32*(t(j))^5;
    Uapprox(n,j) = (1+2*t(j))^5;
end

for i=1:n
    Uapprox(i,2)=10*k*(x(i))^4+Uapprox(i,1);
end

for i = 2:n-1
    for j = 2:m-1
        Uapprox(i, j+1) = (lambda)^2*(Uapprox(i+1,j)-2*Uapprox(i,j)+Uapprox(i-1,j))+2*Uapprox(i,j)-Uapprox(i,j-1);
    end
end

%Here, we plot a graph of the approximate solution with respect to
%both time and space. 

mesh(x,t,Uapprox')
xlabel('$x$', 'interpreter', 'latex')
ylabel('$t$', 'interpreter', 'latex')
zlabel('$u_{approx}$', 'interpreter', 'latex')
title('Approximate Solution with respect to Time and Space')

hold on

