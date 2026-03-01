# scripts/plot_latest.gp
set datafile separator ","
set key left top
set grid

# newest results directory (timestamp format is lexicographically sortable)
latest = system("ls -1 results | sort | tail -n 1")
file = sprintf("results/%s/solution.csv", latest)

# read header to detect optional columns
header = system(sprintf("head -n 1 %s", file))
has_exact = (strstrt(header, "u_exact") > 0)
has_error = (strstrt(header, "error") > 0)

set title sprintf("FEMFramework: %s", latest)
set xlabel "x"

# Plot 1: solution
set ylabel "u"
if (has_exact) {
    plot file using 2:3 with lines title "u_h", \
         file using 2:4 with lines title "u_exact"
} else {
    plot file using 2:3 with lines title "u_h"
}

pause 0.1

# Plot 2: error (if present)
if (has_error) {
    set title sprintf("FEMFramework Error: %s", latest)
    set ylabel "error"
    plot file using 2:5 with lines title "u_h - u_exact"
} else {
    print "No error column found (write_solution_csv without u_exact)."
}

pause -1 "Press Enter to close"