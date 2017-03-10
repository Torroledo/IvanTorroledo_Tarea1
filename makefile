all: time.png

time.png: energy.png
	python timePlot.py
energy.png: energy.dat
	python energyPlot.py
energy.dat: solid.out
	./solid.out 1 > time1.dat
	./solid.out 2 > time2.dat
	./solid.out 4 > time4.dat
solid.out: solid.c
	gfortran solid.c -fopenmp -o solid.out
clean:
	rm *.out
	rm *.dat
	rm *.png
