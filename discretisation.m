function [x,A,b,Mt,map] = discretisation(n,D,mu,L,N,Lbnd,Rbnd)
% Peforms spatial discretisation using finite volume method.

a0 = 1; b0 = 0; aL = 0; bL = 1;

x = linspace(0,L,N); % uniform grid
h = x(2)-x(1); % grid spacing
xw = max(x-h/2,0); % east faces for CVs
xe = min(x+h/2,L); % west faces for CVs
map = @(j,i) (i-1)*N+j; % index for node j of species i

% Source term
f = cell(n,1);
for i = 1:n
    f{i} = zeros(N,1);
end

A = zeros(n*N,n*N);
b = zeros(n*N,1);

% Left boundary node
for i = 1:n
    % species i
    A(map(1,i),map(1,i)) = a0; b(map(1,i)) = Lbnd{i};
end

% Right boundary node
% species 1
A(map(N,1),map(N,1)) = -(D(xw(N))/h + D(xe(1))*aL/bL + mu(1)*h/2);
A(map(N,1),map(N-1,1)) = D(xw(N))/h;
b(map(N,1)) = -D(xe(N))*Rbnd{1}/bL + f{1}(N)*h/2;
% species 2,...,n
for i = 2:n
    A(map(N,i),map(N,i)) = -(D(xw(N))/h + D(xe(1))*aL/bL + mu(i)*h/2);
    A(map(N,i),map(N-1,i)) = D(xw(N))/h;
    A(map(N,i),map(N,i-1)) = mu(i-1)*h/2;
    b(map(N,i)) = -D(xe(N))*Rbnd{i}/bL + f{i}(N)*h/2;
end

% Internal nodes
for j = 2:N-1
    % species 1
    A(map(j,1),map(j-1,1)) = D(xw(j))/h;
    A(map(j,1),map(j,1)) = -(D(xe(j))/h + D(xw(j))/h + mu(1)*h);
    A(map(j,1),map(j+1,1)) = D(xe(j))/h;
    b(map(j,1)) = f{1}(j)*h;
    % species 2
    for i = 2:n
        A(map(j,i),map(j-1,i)) = D(xw(j))/h;
        A(map(j,i),map(j,i)) = -(D(xe(j))/h + D(xw(j))/h + mu(i)*h);
        A(map(j,i),map(j+1,i)) = D(xe(j))/h;
        A(map(j,i),map(j,i-1)) = mu(i-1)*h;
        b(map(j,i)) = f{i}(j)*h;
    end
end

% Divide by CV widths
for i = 1:n
    A(map(1,i),:) = A(map(1,i),:)/(h/2); b(map(1,i)) = b(map(1,i))/(h/2);
    for j = 2:N-1
        A(map(j,i),:) = A(map(j,i),:)/h; b(map(j,i)) = b(map(j,i))/h;
    end
    A(map(N,i),:) = A(map(N,i),:)/(h/2); b(map(N,i)) = b(map(N,i))/(h/2);
end

A = sparse(A); b = sparse(b);

%% Mass matrix for transient problem
Mt = eye(n*N);
for i = 1:n
    Mt(map(1,i),map(1,i)) = 0;
end
Mt = sparse(Mt);
% To match with Eq (14):
A = -A; b = -b;
