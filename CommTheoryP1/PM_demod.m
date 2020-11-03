function result = PM_demod(m,s, f_cutoff, fs)
    Ac = 1.05;
    
    y = Ac*m.*s;
    y_low = lowpass(y,f_cutoff, fs);
    result = downsample (y_low, 40);
end