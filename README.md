# Introduction 
This project examines numerical solutions for the 2-dimensional heat equation using explicit and implicit finite difference methods. explict_dirichlet.m and implicit_dirichlet.m are MATLAB code for explicit and implicit method with homogenous dirichlet condition. 

# Requirement
- MATLAB (version R2018b or later is recommended)
- Open explict_dirichlet.m and implicit_dirichlet.m in MATLAB to run the simulation.

# Parameters
```matlab
L = 1; % Length of the domain (square domain)
N = 51; % Number of grid points in both x and y directions
alpha = 0.1; % Thermal diffusivity
dh = L / (N-1); % Grid spacing (dx = dy = dh)
dt = 0.001; % Time step
gamma = dt*alpha/dh^2;
n_steps = 500; % Number of time steps
```

# Discretization
```matlab
x = linspace(0, L, N);
y = linspace(0, L, N);
[X, Y] = meshgrid(x, y);
```

# Initial conditions
```matlab
initial_condition_func = @(X, Y) 100*sin(pi * X) .* sin(pi * Y); % Initial condition
exact_solution = @(X, Y, t) 100*sin(pi * X) .* sin(pi * Y)*exp(-2*alpha*pi^2*t); % Exact solution
T = initial_condition_func(X, Y);
```

# Boundary condition
```matlab
% Boundary condition
T(1,:) = 0;
T(N,:) = 0;
T(:,1) = 0;
T(:,N) = 0;
```

# Explicit method
```matlab
for n = 1:n_steps
    T_new = T;
    for i = 2:N-1
        for j = 2:N-1
            T_new(i,j) =  (1 - 4*gamma)*T(i,j) + gamma*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1));
        end
    end
    T = T_new;

    % Calculate the exact solution
    T_exact = exact_solution(X, Y, n * dt);

    % Calculate the percentage difference
    percentage_diff = abs((T(2:end-1, 2:end-1) - T_exact(2:end-1, 2:end-1)) ./ T_exact(2:end-1, 2:end-1)) * 100;
  
    error(n) = max(percentage_diff, [], 'all');
end
```

# Implicit method
## Construct the coefficient matrix
```matlab
A = speye(N^2);

for i = 2:N-1
    for j = 2:N-1
        l = (i-1)*N+j;
        A(l, l) = 1 + 4*gamma;
        A(l, l-1) = -gamma;
        A(l, l+1) = -gamma;
        A(l, l-N) = -gamma;
        A(l, l+N) = -gamma;
    end
end
```
## Solving system of linear equation
```matlab
for n = 1:n_steps
    % Reshape T to a column vector
    T_vec = T(:);

    % Solve the linear system
    T_vec_new = A \ T_vec;

    % Reshape the solution back to a matrix
    T = reshape(T_vec_new, N, N);

    % Calculate the exact solution
    T_exact = exact_solution(X, Y, n*dt);

    % Calculate the percentage difference
    percentage_diff = abs((T(2:end-1, 2:end-1) - T_exact(2:end-1, 2:end-1)) ./ T_exact(2:end-1, 2:end-1))*100;
    error(n) = mean(percentage_diff, 'all');
end
```

# visualization
```matlab
imagesc(x, y, T);
title(['Temperature Distribution at Time t = ', num2str(n*dt)]);
xlabel('X');
ylabel('Y');
colorbar;
clim([0, 100]); % Keep the temperature scale constant
axis equal tight;
drawnow;
```

# Plotting error
```matlab
plot(1:n_steps, error);
xlabel('Time Steps');
ylabel('Percentage');
title('Error over Time');
```


