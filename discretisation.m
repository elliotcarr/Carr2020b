function [x,A,b,Mt,map] = discretisation(n,D,mu,L,N,Lbnd)
% Peforms spatial discretisation using finite volume method.

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

%% Left boundary node
% Dirichlet condition
for i = 1:n
    A(map(1,i),map(1,i)) = 1; b(map(1,i)) = Lbnd{i};
end

%% Right boundary node
for i = 1:n
    A(map(N,i),map(N,i)) = -(D(xw(N))/h - mu(i,i)*h/2);
    A(map(N,i),map(N-1,i)) = D(xw(N))/h;
    % reactions
    for k = [1:i-1,i+1:n]
        A(map(N,i),map(N,k)) = mu(i,k)*h/2;
    end
    b(map(N,i)) = -D(xe(N))*0 + f{i}(N)*h/2;
end

%% Internal nodes
for j = 2:N-1
    for i = 1:n
        A(map(j,i),map(j-1,i)) = D(xw(j))/h;
        A(map(j,i),map(j,i)) = -(D(xe(j))/h + D(xw(j))/h - mu(i,i)*h);
        % reactions
        for k = [1:i-1,i+1:n]
            A(map(j,i),map(j,k)) = mu(i,k)*h;
        end
        A(map(j,i),map(j+1,i)) = D(xe(j))/h;
        b(map(j,i)) = f{i}(j)*h;
    end
end

%% Divide by CV widths
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
