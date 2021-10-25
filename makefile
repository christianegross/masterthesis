#Variablen
LIBS:= -lgsl -lgslcblas -lm -fopenmp -lrt

u1.exe: u1.o
	gcc -std=c99 -Wall -pedantic -o $@ $^ $(LIBS)

	
#erstellt aus allen .c dateien eine .o datei	
%.o: %.c
	gcc -std=c99 -Wall -pedantic -fopenmp -I /usr/include/ $^ -c 
	
#.PHONY: clean

#clean:

#$^: fügt alle abhängigkeiten ein
#%: sucht alle Dateien mit gegebener Endung(wildcard)
#$@: Fügt target ein
