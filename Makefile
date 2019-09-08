all: Main clean

Main: Main.o pk_sk_param.o isogeny.o montgomery.o curve_point.o fp.o
	gcc  -g -o Main *.o -lm -lgmp

%.o: %.c
	gcc -c -Wall -g $< -o $@

clean: 
	rm -f *.o
	echo Clean done
