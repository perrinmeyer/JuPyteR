using FFTW
using DSP


function naive_convolution(x,h)

    L = length(x)
    P = length(h)

    y = zeros(L+P-1)
    xp = zeros(L+P-1)
    hp = zeros(L+P-1)
    
    xp[1:L] = x
    xp[L+1:L+P-1] .= 0

    hp[1:P] = h
    hp[P+1:L+P-1] .= 0

    X = fft(xp)
    H = fft(hp)
    Y = X .* H
    y = ifft(Y)

    return(real(y))
    

end

N = 50000

x = rand(N)

y_OVS = zeros(N) 

P = 5

h = rand(P)

L = 17

chunklength = L - P + 1
chunk = zeros(chunklength)

xr = zeros(L) 

yrp = zeros(L+P-1)
yr = zeros(L+P-1)

yconv = conv(x,h)

ynaiveconv = naive_convolution(x,h)

R = Int(floor(N / chunklength)) - 2 

for r=1:R
    for n=1:(L)
        xr[n] = x[n + r * (L - P + 1) - P + 1]
    end
    global yrp = naive_convolution(xr,h)
    for n=(P):(L)
        global yr[n] = yrp[n]
    end
    global chunk = yr[P:L]
    LL = length(chunk)
    global y_OVS[1 + ((r-1) * LL):LL + ((r-1) * LL)] = chunk
    
    
end


tmpstart = L - P +2
ytmp = ynaiveconv[tmpstart:tmpstart+(R * chunklength)-1]

yOVStmp = y_OVS[1:R*chunklength]

ytmp - yOVStmp

maximum(ytmp - yOVStmp)











