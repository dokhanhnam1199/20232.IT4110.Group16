% Parameters
L = 1; % Length of the domain (square domain)
N = 51; % Number of grid points in both x and y directions
alpha = 0.1; % Thermal diffusivity
dh = L / (N-1); % Grid spacing (dx = dy = dh)
dt = 0.001; % Time step
gamma = dt*alpha/dh^2;
n_steps = 500; % Number of time steps

% Discretized spatial domain
x = linspace(0, L, N);
y = linspace(0, L, N);
[X, Y] = meshgrid(x, y);

% Initial temperature distribution
initial_condition_func = @(X, Y) 100*sin(pi * X) .* sin(pi * Y); % Initial condition
exact_solution = @(X, Y, t) 100*sin(pi * X) .* sin(pi * Y)*exp(-2*alpha*pi^2*t); % Exact solution
T = initial_condition_func(X, Y);

% Boundary condition
T(1,:) = 0;
T(N,:) = 0;
T(:,1) = 0;
T(:,N) = 0;

% Error
error = zeros(1, n_steps);

% Start timing
tic;

% Construct the coefficient matrix
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

% Create a figure for the plot
figure;

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

    % Plot the temperature distribution
    imagesc(x, y, T);
    title(['Temperature Distribution at Time t = ', num2str(n*dt)]);
    xlabel('X');
    ylabel('Y');
    colorbar;
    clim([0, 100]); % Keep the temperature scale constant
    axis equal tight;
    drawnow;
end

% Stop timing
time = toc;

% Display the running time
fprintf('Running time for Implicit method: %.6f seconds\n', time);

% Plot error difference over time
figure(2);
plot(1:n_steps, error);
xlabel('Time Steps');
ylabel('Percentage');
title('Error over Time');