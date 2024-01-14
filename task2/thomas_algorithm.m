function U = thomas_algorithm(A, B, C, G)
    
    N = size(B, 2);
    S = zeros(N - 1, 1);
    T = zeros(N, 1);
    U = zeros(N, 1);
    S(1) = -C(1) / B(1);
    T(1) = G(1) / B(1);

    for i = 2:1:(N - 1)
        S(i) = -C(i) / (B(i) + A(i) * S(i - 1));
    end
    for i = 2:1:N
        T(i) = (G(i) - A(i) * T(i - 1)) / (B(i) + A(i) * S(i - 1));
    end
    U(N) = T(N);

    for i = (N - 1):-1:1
        U(i) = S(i) * U(i + 1) + T(i);
    end
end