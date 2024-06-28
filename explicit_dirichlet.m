% Parameters
L = 1; % Length of the domain (square domain)
N = 51; % Number of grid points in both x and y directions
alpha = 0.1; % Thermal diffusivity
gamma = 0.25;
dh = L / (N-1); % Grid spacing (dx = dy = dh)
dt = gamma * dh^2 / alpha; % Time step (stability condition for explicit method)
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

% time-stepping loop
figure;

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
fprintf('Running time for Explicit method: %.6f seconds\n', time);

% Plot error difference over time
figure(2);
plot(1:n_steps, error);
xlabel('Time Steps');
ylabel('Percentage');
title('Error over Time');