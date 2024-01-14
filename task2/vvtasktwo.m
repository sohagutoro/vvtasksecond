format shortE

dt = 0.001;
dx = 0.05;
X_max = 1;
T_max = 5;

X = 0:dx:X_max;
T = 0:dt:T_max;
U1 = zeros(size(T, 2), size(X, 2));
U2 = U1;
U3 = U1;

f = @(x, t) exp(-t) .* (x.^2 - x + 2);
g = @(x) cos(2 * x) + (1 - x) .* x;
a = @(t) exp(-4 * t);
b = @(t) exp(-4 * t) * cos(2);
u = @(x, t) exp(-4 * t) .* cos(2 * x) + exp(-t) .* (1 - x) .* x;

U1(1, :) = g(X);
U2(1, :) = U1(1, :);
U3(1, :) = U1(1, :);

[X_, T_] = meshgrid(X, T);
U = u(X_, T_)

for i = 2:1:size(T, 2)
    U1(i, 2:end-1) = U1(i-1, 2:end-1) + dt * (f(X(2:end-1), T(i-1)) + L(U1(i-1, :), dx));
    U1(i, 1) = a(T(i));
    U1(i, end) = b(T(i));
end

U1

for i = 2:1:size(T, 2)
    % Расчет трех диагоналей и столбца свободных членов
    A = [0, 1 / dx ^ 2 * ones(1, size(U2, 2)-2), 0];
    B = [1, -(2 / dx ^ 2 + 1 / dt) * ones(1, size(U2, 2)-2), 1];
    G = [a(T(i)), -1 / dt * U2(i-1, 2:end-1) - f(X(2:end-1), T(i)), b(T(i))];
    % Решение системы методом прогонки
    U2(i, 1:end) = thomas_algorithm(A, B, A, G);
end

U2

for i = 2:1:size(T, 2)
    % Расчет трех диагоналей и столбца свободных членов
    A = [0, 0.5 / dx ^ 2 * ones(1, size(U3, 2)-2), 0];
    B = [1, -(1 / dx ^ 2 + 1 / dt) * ones(1, size(U3, 2)-2), 1];
    G = [a(T(i)), -1 / dt * U3(i-1, 2:end-1) - f(X(2:end-1), T(i) - dt / 2) - ...
        0.5 * L(U3(i-1, :), dx), b(T(i))];
    % Решение системы методом прогонки
    U3(i, 1:end) = thomas_algorithm(A, B, A, G);
end

U3

max(abs(U1 - U), [], 'all')

max(abs(U2 - U), [], 'all')

max(abs(U3 - U), [], 'all')

figure
hold on

% Построение трехмерного точечного графика для U1
figure;
scatter3(X_, T_, U1, 20, 'r', 'filled');
title('График точек Явной схемы');
xlabel('x');
ylabel('t');
zlabel('U1');
grid on;

% Построение трехмерного точечного графика для U2
figure;
scatter3(X_, T_, U2, 20, 'g', 'filled');
title('График точек Неявной схемы');
xlabel('x');
ylabel('t');
zlabel('U2');
grid on;

% Построение трехмерного точечного графика для U3
figure;
scatter3(X_, T_, U3, 20, 'b', 'filled');
title('График точек Кранка-Никольсона');
xlabel('x');
ylabel('t');
zlabel('U3');
grid on;