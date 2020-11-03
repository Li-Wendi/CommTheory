function result = SSB_demod(m, c, f_cutoff, fs)
    m= m.*c;
    result = lowpass(m, f_cutoff, fs);
    result = highpass(result, 2 ,fs);
    result = downsample(result, 40);
end