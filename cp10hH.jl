
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





