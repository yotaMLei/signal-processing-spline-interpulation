function x = interpCubic(c,t)
% interpolation using spline of order 3

% filter using the coefficients we calculated.
a_C = [1 2-sqrt(3)];% coefficients for the causal filter.
a_NC = [(2+sqrt(3))/6 1/6];% coefficients for the non-causal filter.

d = filter(1, a_C, c);% applying the causal filter.
d = filter(1, a_NC, wrev(d));% applying the non-causal filter on the reversed signal.
d = wrev(d);% reverse back.
d = d(2:end);

x = SplineExpansion(d,t,3);% calculate the interpulation with spline of order 3.

end

