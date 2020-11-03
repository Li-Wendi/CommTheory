function result = FM_mod(m, t, fc, fs, Kf)
    Ac = 11/5;
    result = Ac*cos(2*pi*fc*t + 2*pi*Kf*cumsum(m)/fs);
end

