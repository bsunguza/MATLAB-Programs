function  SunguzaFDMPoissonDirichelet

%The following code will numerically solve a Poisson PDE Equation within a unit square
%whose bottom left corner is on the origin. Using latex notation, the
%equation is $-\nabla u = 2\pi^2 \sin(\pi x)\sin(\pi y)$ where we have 
%homogenous Dirichelet boundary conditions. 

%In this code, we will consider four step sizes that will be generated 
%by the number of steps N=10,20,40, and 80. After finding a numerical 
%approximation of the solution of the above PDE, we will consider the 
%maximum error of the approximation as a function of h, the step size. 

N = 5; 
maxEh = zeros(4,1);
H = zeros(4,1);

for r = 1:4

N = 2*N; 
h = 1/N;

H(r) = h;

x = 0:h:1;
y = 0:h:1;

%The code that follows is based on an email from Kening Wang 
%that detailed how to make a block matrix of the kind found in this problem. 
%Note that we are currently focusing on the unknowns of the problem, 
%that is, the interior points. Hence, we are going to approximate the 
%solution on a N-2 x N-2 grid.

nBlocks = N-2; % number of blocks

Q = diag(-4*ones(N-2,1))+diag(ones(N-3,1), -1)+diag(ones(N-3,1), 1); % diagonal block
Q(1, N-2) =1;
Q(N-2, 1)=1;
Up = eye(N-2); % upper diagonal block
Low = eye(N-2); % lower diagonal block
Top = zeros(N-2); Top(1,N-2) =1; %Upper corner
Bottom = zeros(N-2); Bottom(N-2,1)=1; %Lower corner
I = eye(nBlocks); % block identity
E = diag(ones(nBlocks-1,1),1); % upper diagonal (1's)
F = diag(ones(nBlocks-1,1),-1); % lower diagonal (1's)
A = kron(I, Q) + kron(E, Up) + kron(F, Low) + kron(Top,I)+kron(Bottom,I);

%After noticing that A was illconditioned, I decided to change one of the
%entries a viable approximation. This change is shown on line 52. The maximum 
%deviation per h of the resulting approximations are captured within the 
%error Eh that is calculated starting on line 131. 

A(1,1) = 1;

%The next section is dedicated to obtain the vector b in the matrix-vector
%equation Au = b. Since there are two indices involved, I decided to use
%counters external to the 'for loop' below in order to increment the
%indices as I saw fit. 

b = zeros((N-3)^2,1);
j = 2;
i = 2;
multiple = 1;

for k = 1:(N-2)^2

    b(k) = -2*((pi)^2)*(h^2)*sin(pi*x(i))*sin(pi*y(j));
    i = i+1;

    if k == multiple*(N-2) 
        i = 2;
        j = j+1;
        multiple = multiple+1;
    end
end

%Next, we determine the approximation Uapprox. Note that Psi is the
%approximation of the interior points. The 'for loop' below is used to find
%a complete approximation for the interior points as well as all of the
%boundaries. Again, counters were used to assist in the assignment of
%values. 

Psi = A\b;

Uapprox = zeros(N*N,1);
multiple = 1;
counter = 2;

for i = N:(N-1)*N

    if i == multiple*N
        Uapprox(i) = 0;
    end

    if i == multiple*N+1
        multiple = multiple+1;
        counter = counter-1;
        Uapprox(i) = 0;
    elseif i < multiple*N
        Uapprox(i) = Psi(counter);
        counter = counter+1;
    end
end

%Now, we determine Uexact, using the techniques that were used above. Note
%that the classical solution is known for this problem and in latex notation
%it is given by $u(x,y) = \sin(\pi x)\sin(\pi y)$ on the entire unit square
%described in the first comment. 

Uexact = zeros(N^2,1);

j = 1;
i = 1;
multiple = 1;

for k = 1:(N)^2

    Uexact(k) = sin(pi*x(i))*sin(pi*y(j));
    i = i+1;

    if k == multiple*N 
        i=1;
        j = j+1;
        multiple = multiple+1;
    end
end

Uexact;

%With Uapprox and Uexact, we can now determine the absolute errors as well
%as the maximum absolute error for all points in the grid of interest. 

Eh = zeros(N*N,1); 

for i = 1:N*N
    Eh(i) = abs(Uapprox(i) - Uexact(i));
end

Ehmax = Eh(1);

for i = 1:N*N
    if Ehmax<Eh(i)
        Ehmax = Eh(i);
    end
end

maxEh(r) = Ehmax;

end

%Finally, we use a log-log plot to look at the maximum error as a function
%of h. 

loglog(H,maxEh)
xlabel('$h$', 'interpreter', 'latex')
ylabel('$E_h$', 'interpreter', 'latex');
title('log-log plot of error $E_h$ and step size $h$', 'interpreter', 'latex')

end