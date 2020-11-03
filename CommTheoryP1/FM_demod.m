function [result] = FM_demod(m, fs, f_cutoff, w)
w = w*2*pi;
y = real(ifft(((fft(m)/fs)).* (j * w)));

y(y < 0) = 0;  

y_low = lowpass(y,f_cutoff,fs);

result = decimate(y_low,40);

end