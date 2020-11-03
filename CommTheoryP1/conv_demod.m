function result = conv_demod(m, fs, f_cutoff)
    m(m<0) = 0;
    y = lowpass(m, f_cutoff, fs);
    result = downsample(y, 40);
end