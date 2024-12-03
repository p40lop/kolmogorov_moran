reset
set term pngcairo size 600,400 font "FreeSans,12"
set encoding utf8
# set xrange [0:7]
set xlabel "y"
set xrange [0:15]
set ylabel "u_x(y)" rotate by 0
set format y "%+2.1t \U+00B7 10^{%T}"
set ytics add ("0" 0)
nu(tau)=(tau-0.5)/(3)
u_x_teo(x,tau)=1e-7/nu(tau)*(16/6.28)**2*sin(6.28*x/16)


#WIN
if (strstrt(GPVAL_SYSNAME, "Windows") == 1) {
    command = "powershell New-Item -ItemType Directory -Force -Path "
}
else {
    command = "mkdir -p "
}

system(command."'figures'")
system(command."'figures/test'")

taus= system("ls data")
do for [tau in taus]{
    files=system("ls data/".tau."/time/")
    system(command."'figures/test/".tau."'")
    system(command."'figures/test/".tau."/time'")

    # A plot for each saved time to check simmetry

    do for [file in files]{
        set title gprintf("tau = %4.2f ",real(tau[5:])).gprintf("t = %+2.1t \U+00B7 10^{%T}",real(file[:6]))
        set out "figures/test/".tau."/time/".tau."_t_".file[:6].".png"
            plot "data/".tau."/time/".file u 2:($1==7?$4:1/0)  w p pt 7 t "x=7",\
                 ""                   u 2:($1==4?$4:1/0)  w p pt 7 t "x=4" ,\
                 ""                   u 2:($1==1?$4:1/0)  w p pt 7 t "x=1" ,\
                 u_x_teo(x,real(tau[5:])) t "u_x^{teo}" lt -1 lw 2
        unset out

    }
}

set term pngcairo size 600,600 font "FreeSans,12"
# A plot to show time evolution
set origin 0, 0.37
set size 1,0.65
do for [tau in taus]{
    files=system("ls data/".tau."/time/")
    set key at screen 0.5, screen 0.2 center title "t" vertical
    set title gprintf("tau = %4.2f",real(tau[5:]))
    set out "figures/test/".tau."/".tau."_t_evo.png"
    plot for [file in files] "data/".tau."/time/".file u 2:($1==3?$4:1/0)  w p pt 7 t gprintf("%+2.1t \U+00B7 10^{%T}",real(file[:6])),\
                 u_x_teo(x,real(tau[5:])) t "u_x^{teo}" lt -1 lw 2
    unset out
    set key inside
}

# MORAN
set term pngcairo size 600,400 font "FreeSans,12"
set origin 0,0
set size 1,1
unset xrange
unset yrange
unset key
set log x
set xlabel "t"
set format x "10^{%T}"
set format y "%+3.2t \U+00B7 10^{%T}"
do for [tau in taus]{
    set title gprintf("tau = %4.2f",real(tau[5:]))
    set out "figures/test/".tau."/".tau."_I_t.png"
    plot "data/".tau."/macro.dat" u 1:2 w lp pt 7 t ""
    unset out
}

set xrange [0:3]
set yrange [0:15]
set size ratio 1
set format xy "%2.0f"
set xlabel "x"
set ylabel "y"
set cblabel "I_{loc}" rotate by 0 offset screen 0.01
do for [tau in taus]{

    system(command."'figures/test/".tau."/moran'")


    files=system("ls data/".tau."/time/")
    # A plot for each saved time
    do for [file in files]{
        set title gprintf("tau = %4.2f ",real(tau[5:])).gprintf("t = %+2.1t \U+00B7 10^{%T}",real(file[:6]))
        set out "figures/test/".tau."/moran/".tau."_moran_t_".file[:6].".png"
        plot "data/".tau."/time/".file u 1:2:6 w image t ""
        unset out


        }
    }
