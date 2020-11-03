function result = SSB_mod (mt, s, c,f);
    M = fft(mt);
    Ac = 1
    H_m = Hilbert_Function (M, f);
    result = Ac*(mt.*c + H_m.*s) + c; 
end