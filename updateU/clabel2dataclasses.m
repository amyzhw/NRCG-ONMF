function dataclasses = clabel2dataclasses(clabel, nc)
n = length(clabel);
dataclasses = full(sparse((1:n)', clabel(:), ones(n,1), n, nc));
