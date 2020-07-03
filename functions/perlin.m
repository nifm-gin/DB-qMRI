function s = perlin(m)
  s = zeros(m);    % output image
  w = m;           % width of current layer
  i = 0;           % iterations
  while w > 3
    i = i + 1;
    d = interp2(randn(w), i-1, "spline");
    s = s + i * d(1:m, 1:m);
    w = w - ceil(w/2 - 1);
  end
end