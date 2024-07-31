function out = getavgLineout(tt,wl,data,tq,Lq,tstp,N)
    out = 0;
    for ii = 0:(2 * N)
        out = out + interp2(tt,wl,data,tq-(N-ii)*tstp,Lq); 
    end
    out = out / (2 * N);
end