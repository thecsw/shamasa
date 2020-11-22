sigma = 10;
beta = 8/3;
rho = 28;

f = @(t,a) [
    -sigma*a(1) + sigma*a(2);
    rho*a(1) - a(2) - a(1)*a(3);
    -beta*a(3) + a(1)*a(2)];

[t,a] = ode45(f, 0:0.01:100, [1 1 1]);

% Plot the ideal Lorenz system
figure('name', 'Lorenz');
plot3(a(:,1),a(:,2),a(:,3), '-', 'Color', '#B28DFF');
title('Lorenz 63, timestep: h=0.01, range: 0->100');
xlabel('X');
ylabel('Y');
zlabel('Z');

% Let's add some noise that is centered around zero
% Adjust this value to manipulate Sigma for mvnrnd
noise = 0.05;
obs = a + mvnrnd(zeros(size(a)), noise*ones(1, width(a)));
figure('name', 'Noisy Lorenz');
plot3(obs(:,1),obs(:,2),obs(:,3), '-', 'Color', '#FFABAB');
title('Lorenz 63, timestep: h=0.01, range: 0->100');
xlabel('X');
ylabel('Y');
zlabel('Z');

% Now, let's try to build the Runge Kutta

