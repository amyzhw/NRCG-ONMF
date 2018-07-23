function i = ip(x,h1,h2)
%% IP	computes the inner product of h1,h2 which are tangents at x


i = h1.Up(:)'*h2.Up(:) + h1.Vp(:)'*h2.Vp(:);
