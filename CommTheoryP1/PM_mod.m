function[result] = PM_mod(m,T, fc,kf)
    Ac =1.05;
    N = length(m);
    t = linspace(0,T,N);
    result = Ac*cos(2*pi*fc*t + kf*m);
end