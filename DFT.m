% Density Functional Theory (DFT) project

% Define constants
h_bar = 1.0545718e-34; % Planck's constant / (2 * pi) in J*s
mass_e = 9.10938356e-31; % Electron mass in kg
e = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = 8.854187817e-12; % Vacuum permittivity in F/m
pi = 3.14159265359; 

% Define parameters
N = 100; % Number of grid points
L = 1e-9; % Length of the box in metres
dx = L / N; % Grid spacing
x = linspace(0, L, N); % Grid points

% Define the potential energy function
V_ext = zeros(1, N); % External potential (for simplicity, set to zero)

% Set up the electron density (initial guess)
n = ones(1, N) * 2; % Electron density (two electrons per site for simplicity)

% Initialise the loop
max_iterations = 100; % Maximum number of iterations
tolerance = 1e-6; % Tolerance for convergence
converged = false;
iteration = 0;

while ~converged && iteration < max_iterations
    % Solve the Kohn-Sham equations
    % Using a simple finite difference approximation
    
    % Calculate the effective potential
    V_eff = V_ext + (e^2 / (4 * pi * epsilon_0)) * sum(n .* (1./sqrt((x - x').^2 + dx^2) - eye(N)), 2) * dx;
    
    % Solve the Schrödinger equation for each electron
    [eigenvecs, ~] = eig(-h_bar^2 / (2 * mass_e) * second_derivative(N, dx) + diag(V_eff), 'vector');
    
    % Update the electron density
    n_new = 2 * sum(abs(eigenvecs).^2, 2)';
    
    % Check for convergence
    if max(abs(n_new - n)) < tolerance
        converged = true;
    end
    
    % Update electron density
    n = n_new;
    
    % Increase iteration count
    iteration = iteration + 1;
end

if converged
    disp('Converged');
else
    disp('Did not converge within maximum iterations');
end


% Calculate total energy
kinetic_energy = sum(sum(eigenvals .* abs(eigenvecs).^2)) * (-h_bar^2 / (2 * mass_e));
total_energy = kinetic_energy + sum(n .* V_eff) * dx;

disp(['Total Energy: ', num2str(total_energy), ' J']);

% Plots
figure;
subplot(2,1,1);
plot(x, n, 'b-', 'LineWidth', 1.5);
xlabel('Position (m)');
ylabel('Electron Density');
title('Electron Density');

subplot(2,1,2);
plot(x, V_eff, 'r-', 'LineWidth', 1.5);
xlabel('Position (m)');
ylabel('Effective Potential (J)');
title('Effective Potential');

% Function to calculate second derivative using finite differences
function D2 = second_derivative(N, dx)
    D2 = (-2 * diag(ones(1, N)) + diag(ones(1, N-1), 1) + diag(ones(1, N-1), -1)) / (dx^2);
end

