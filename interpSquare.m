function x = interpSquare(c,t)
% interpolation using spline of order 2

% filter using the coefficients we calculated.
a_C = [1 3-2*sqrt(2)];% coefficients for the causal filter.
a_NC = [(3+2*sqrt(2))/8 1/8];% coefficients for the non-causal filter.

d = filter(1, a_C, c);% applying the causal filter.
d = filter(1, a_NC, wrev(d));% applying the non-causal filter on the reversed signal.
d = wrev(d);% reverse back.
d = d(2:end);

x = SplineExpansion(d,t,2);% calculate the interpulation with spline of order 2.

end

