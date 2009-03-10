function y = checkergauss(N)

hN = ceil(N/2);
y = zeros(N);
for i = 1:N
    for j = 1:N
        g = exp(-(((i-hN)/hN)^2 + (((j-hN)/hN)^2))*4);
        if xor(j>i,j>N-i)
            y(i,j) = -g;
        else
            y(i,j) = g;
        end
    end
end