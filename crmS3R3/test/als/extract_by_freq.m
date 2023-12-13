function sig = extract_by_freq(sig, nmin_f, nmax_f)

gis = fft(sig);

if isfinite(nmin_f)
   gis(1:nmin_f-1) = 0;
   gis(end-nmin_f+2:end) = 0;
end

if isfinite(nmax_f)
   gis(nmax_f+1:end-nmax_f-1) = 0;
end

sig = real(ifft(gis));

