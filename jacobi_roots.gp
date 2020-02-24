gnuplot

set xlabel "Re(x)"
set ylabel "Im(x)"

set title "Jacobi polynomail roots P_{50}^{(2,-42.5)}(x)"
plot 'jacobi_roots.dat' index 0

set title "Jacobi polynomail roots P_{50}^{(2,-52)}(x)"
plot 'jacobi_roots.dat' index 1

set title "Jacobi polynomail roots P_{50}^{(2,-63.5)}(x)"
plot 'jacobi_roots.dat' index 2


set title "Jacobi polynomail roots P_{50}^{(-42.5,2)}(x)"
plot 'jacobi_roots.dat' index 3

set title "Jacobi polynomail roots P_{50}^{(-52,2)}(x)"
plot 'jacobi_roots.dat' index 4

set title "Jacobi polynomail roots P_{50}^{(-63.5,2)}(x)"
plot 'jacobi_roots.dat' index 5
