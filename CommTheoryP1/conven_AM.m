function [result] = conven_AM( m, fc, t,a) 
c = cos(2*pi*fc*t);
Ac = 6/5;
result = Ac/2*(1 + a*m).* c;
%conv_AM = fft(conv_am);
%conv_AM_s = fftshift(conv_AM/fs);
end

%sounds worse because rectifiler, it is not ideal rectifier