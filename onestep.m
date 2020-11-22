% ic is the initial point
% h is the step size
function y = onestep(f, h, ic)
    k1 = f(0, ic)';
    k2 = f(0, ic+0.5*h* k1)';
    k3 = f(0, ic+0.5*h*k2)';
    k4 = f(0, ic+h*k3)';
    y = ic+((1/6)*h*(k1+2*k2+2*k3+k4));
end