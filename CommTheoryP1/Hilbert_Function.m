function output = Hilbert_Function(input, freq)
    inside = (-j*sign(freq));
    output = real(ifft(input.*(inside)));
end