
using PyPlot
using DSP 
using FFTW 

function naive_filter(b,x)

    N = length(x)
    M = length(b)
    
    y = zeros(N)
    
    for n=M:N

        for k=0:(M-1)
            y[n] = y[n] + (b[k+1] * x[n-k])
        end
    end
    
    return y 

    
end


function naive_convolution(b,x)

    L = length(x)
    P = length(b)

    y = zeros(L+P-1)
    xp = zeros(L+P-1)
    hp = zeros(L+P-1)
    
    
    xp[1:L] = x
    xp[L+1:L+P-1] .= 0

    hp[1:P] = b
    hp[P+1:L+P-1] .= 0

    X = fft(xp)
    H = fft(hp)
    Y = X .* H
    y = ifft(Y)

    return(real(y))
    

end



f_samp = 48e3

dt = 1 / f_samp


N = 1000

n = collect(0:(N-1))

x1 = sin.(2 * pi * 20 * n ./ f_samp)

x2 = sin.(2 * pi * 20000 * n ./ f_samp)

x = x1 + x2

b = ones(10) / 10
a = zeros(5)
a[1] = 1

y = filt(b,a,x)

yc = conv(b,x)

ynf = naive_filter(b,x)

ync = naive_convolution(b,x)

R = zeros(20,4)


R[:,1] = y[21:40]
R[:,2] = yc[21:40]
R[:,3] = ynf[21:40]
R[:,4] = ync[21:40]









figure(1)
## plot(x[1:1000],"-r")
## plot(y[1:1000],"-g")
plot(x1,"-r")
plot(x2,"-g")
plot(x,"-b")
plot(y,"-k")




    




