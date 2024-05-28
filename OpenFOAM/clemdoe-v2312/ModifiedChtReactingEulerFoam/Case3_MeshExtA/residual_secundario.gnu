set logscale y
set title "Residuals"
set ylabel 'Residual'
set xlabel 'Iteration'
#set yrange [1e-7:]
#set xrange [1132000:]
plot "< cat log |grep 'GAMG:  Solving for p_rgh, Initial residual =' | cut -d' ' -f9 | tr -d =,= " title 'Residual p_rgh' with lines
pause 5
reread
