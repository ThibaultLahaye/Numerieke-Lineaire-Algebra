function a = approxCheby(f,n)
    % Determining the real series f(z0), f(z1), ... , f(zn)
    z = zeros(1,n+1);
    v = zeros(1,n+1);
    for l = 1:n+1
        z(1,l) = cos(pi*(l-1)/n);
        v(1,l) = f(z(l));
    end

    % Calculating factors ak
    w = fliplr(v(2:end-1));
    v_even = [v w];
    V = fft(v_even)/n;
    V = real(V(1:n+1));
    a = [V(1)/2 V(2:end)];
end