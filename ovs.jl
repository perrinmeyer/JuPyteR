
include("simplot.jl")
using PyPlot
using DSP 
using FFTW 
using WAV


function cp10hH(centerHz,octave3dB,gaindB,f_samp,P)
    
    f_nyq = f_samp / 2;
    T = 1 ./ f_samp;

    if gaindB >= 0
        boostcutflag = 1
    else
        boostcutflag = 0
    end
    
    gaindB = abs(gaindB) 

    ## % Pf is prewarped frequency (from Matlab "doc bilinear")
    f_p = centerHz;
    Pf = (2 .* pi .* f_p) / ( tan(pi .* (f_p / f_samp)));
    G = 10^(gaindB / 20);
    Q = 1.43 / octave3dB;

    W = 2 * pi * centerHz;

    bzmath = zeros(3,1)
    azmath = zeros(3,1) 

    if boostcutflag == 0
        ## CUT equation
        ## from Mathematica 
        bzmath[1] = Pf.^2 .* Q .+ Pf .* W .+ Q .* W.^2;
        bzmath[2] = -2 .* Pf.^2 .* Q .+ 2 .* Q .* W.^2;
        bzmath[3] = Pf.^2 .*Q .- Pf .* W .+ Q .* W.^2;

        azmath[1] = Pf.^2 .* Q .+ G .* Pf .* W + Q .* W.^2;
        azmath[2] =  -2 .* Pf.^2 .* Q .+ 2 .* Q .* W.^2;
        azmath[3] =  Pf.^2 .*Q - G .* Pf .* W .+ Q .* W.^2;

        bzmath = bzmath ./ azmath[1];
        azmath = azmath ./ azmath[1];
    end

    
    if boostcutflag == 1 
    
        ## BOOST equation 
        bzmath[1] = Pf.^2 .* Q .+ G .* Pf .* W .+ Q .* W.^2;
        bzmath[2] = -2 .* Pf.^2 .* Q .+ 2 .* Q .* W.^2;
        bzmath[3] = Pf.^2 .* Q .- G .* Pf .* W .+ Q .* W.^2;
  
        azmath[1] = Pf.^2 .* Q .+ Pf .* W .+ Q .* W.^2;
        azmath[2] = -2 .* Pf.^2 .* Q .+ 2 .* Q .* W.^2;
        azmath[3] = Pf.^2 .* Q .- Pf .* W .+ Q .* W.^2;

        bzmath = bzmath ./ azmath[1];
        azmath = azmath ./ azmath[1];

    end


    b0 = bzmath[1]
    b1 = bzmath[2]
    b2 = bzmath[3]

    a0 = azmath[1];
    a1 = azmath[2];
    a2 = azmath[3];

    b = bzmath;
    a = azmath;

    ## calculationg H(e^{j w}) (Fourier transform) 
    n = 0:(P-1);
    n2 = 0:(P/2);
    f2 = (f_samp .* n2) ./ P
    H2 = zeros(div(P,2) + 1,1)

    z = exp.(im .* 2 .* pi .* (n ./ P));
    Htop  = b0 .+ b1 .* z.^(-1) .+ b2 .* z.^(-2);
    Hbot = 1 .+ a1 .* z.^(-1) .+ a2 .* z.^(-2);
    H = Htop ./ Hbot;

    H2 = H[1:div(P,2)+1]

    B = b2/a2;
    c1 = b1 - ((b2/a2)*a1);
    c0 = b0 - (b2/a2);

    alpha = (-a1 + sqrt(complex(a1^2 - 4 * a2)))/2;
    beta = (-a1 - sqrt(complex(a1^2 - 4 * a2)))/2;

    A1 = ((c1 / alpha) + c0) / ( 1 - (beta / alpha));
    
    A2 = ( (c1 / beta) + c0) / ( 1 - (alpha / beta));

    h = zeros(P,1)


    h[n[1]+1] = B .+ A1 .+ A2;
    h[ n[2:end] .+ 1,1] = A1 .* alpha.^( n[2:end]) .+ A2 .* beta.^(n[2:end]);

    ## note equivalence of h versus hifft f
    ##maximum(  abs.(h - real(ifft(H)) ))

    return f2,H2,b,a,h,H
    
end

function naive_convolution(h,x)

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

f2,H2,b,a,h,H= cp10hH(2000,0.5,-28,48e3,2^12);

f_samp = 48e3
dt = 1 / f_samp
N = 50000
n = collect(0:(N-1))
x1 = sin.(2 * pi * 100 * n ./ f_samp)
x2 = sin.(2 * pi * 2000 * n ./ f_samp)
x = x1 + x2;

y = filt(vec(b),vec(a),x);

htrunc = h[1:1024];

a = zeros(length(htrunc))
a[1] = 1.0

ytrunc = filt(vec(htrunc),a,x);

yfft = naive_convolution(vec(htrunc),x);


# Julia Overlap-save algorithm setup 
y_OVS = zeros(N) 
P = length(htrunc)
# our "chunk length" for block convolution is 256
L = P + 256
chunklength = L - P + 1
chunk = zeros(chunklength)

xr = zeros(L) 

yrp = zeros(L+P-1)
yr = zeros(L+P-1)
R = Int(floor(N / chunklength)) - 2 


for r=4:R
    # pick out a "chunk from the input stream x[n]"
    for n=1:(L)
        xr[n] = x[n + r * (L - P + 1) - P + 1]
    end
    # do the FFT (block/chunk) convolution
    global yrp = naive_convolution(xr,htrunc)
    # pick out the valid samples 
    for n=(P):(L)
        global yr[n] = yrp[n]
    end
    global chunk = yr[P:L]
    LL = length(chunk)
    # add the valid samples to the output (filtered) stream y_OVS[]
    global y_OVS[1 + ((r-1) * LL):LL + ((r-1) * LL)] = chunk
end

# our naive block convolution algorithm does not handle startup yet, so pick out valid "middle" 
y_OVS_trunc = y_OVS[772:end-1000]

# pick out the correspoding values from our earlier yfft convolution 
yfft_trunc = yfft[1029:1029 + (length(y_OVS_trunc) -1)]

# compare

maximum(y_OVS_trunc - yfft_trunc)










