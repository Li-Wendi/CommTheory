function result = SSB_mod (mt, s, c,f)
    M = fft(mt);
    Ac = 19/10;
    H_m = Hilbert_Function (M, f);
    result = Ac*(mt.*c - H_m.*s) + c; 
end