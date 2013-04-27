function nbz = nbturn (hh)
% function nbz = nbturn (hh)
% count the number of turns of the mother wavelet defined by hh

[x,p,xx,w] = scalewave(hh,5);
nbz = 0;
for k = 1 : length(w)-2
    if ((w(k) - w(k+1)) * (w(k+1) - w(k+2))) < 0
        nbz = nbz + 1;
    end;
end;