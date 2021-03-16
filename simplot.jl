
function simplot(x,y ;
                 psmaxisrange = [20 , 20e3 , -30 , 30],
                 psmxlabel = "Frequency in Hz",
                 psmylabel = "dB",
                 psmtitle = "SIMplot semilogx and dB",
                 psmlinestyle = "-k",
                 psmticklabelsize = 8,
                 psmaxislabelsize = 10,
                 psmtitlesize = 12,
                 )

#====================================================================
[[ note semicolon ]] 
                       x 
julia> function f(p, q ; r = 4, s = "hello")
  println("p is $p")
  println("q is $q")
  return "r => $r, s => $s"
end
f (generic function with 1 method)
When called, this function expects two arguments, and will also accept
a number and a string, labelled r and s. If you don't supply them, the
default values are used:
====================================================================#
    
    ## matplotlib rcparams (startup...)

    rc("xtick",labelsize=6)
    rc("ytick",labelsize=6)

    psmplot = semilogx(x,20.0 .* log10.(abs.(y)),psmlinestyle,linewidth=3)
    psmaxis = axis(psmaxisrange)
    psmaxes = gca()
    thirdfreqs = ceil.(1000 * (2).^((-17:13)/3));
    psmaxes.set_xticks(thirdfreqs)
    psmaxes.xaxis.set_tick_params(which="minor",length=0,width=0)
    psmaxes.xaxis.set_tick_params(which="major",direction="in",labelsize=psmticklabelsize)
    psmaxes.yaxis.set_tick_params(which="major",direction="in",labelsize=psmticklabelsize)


    ## Xlabels_Dict = matread("Xlabels.mat");
    ## this seems to work 
    Xlabels_Dict =
        Dict("Xlabels" =>
             String[" 20", " 25", " 32", " 40", 
                    "50", " 63", " 80", " 100", " 125", " 160", " 200", " 250", " 315", 
                    "400", " 500", " 630", " 800", "1.0k", "1.3k", "1.6k", "2.0k", "2.5k",
                    "3.2k", "4.0k", "5.0k", "6.3k", "8.0k", " 10k", " 13k", " 16k", 
                    "20k"])

    psmaxes.set_xticklabels(Xlabels_Dict["Xlabels"])
    grid(linestyle=":",color="black",linewidth=0.3)
    xlabel(psmxlabel,fontsize=psmaxislabelsize)
    ylabel(psmylabel,fontsize=psmaxislabelsize)
    title(psmtitle,fontsize=psmtitlesize) 

end
