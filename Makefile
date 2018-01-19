all: Main clean

Main: Main.o pk_sk_param.o isogeny.o montgomery.o curve_point.o fp.o
	gcc  -g -o Main *.o -lm -lgmp

Main.o: Main.c
	gcc -c -Wall -g Main.c

pk_sk_param.o: pk_sk_param.c
	gcc -c -Wall -g pk_sk_param.c

isogeny.o: isogeny.c
	gcc -c -Wall -g isogeny.c 

montgomery.o: montgomery.c
	gcc -c -Wall -g montgomery.c 
 
curve_point.o: curve_point.c
	gcc -c -Wall -g curve_point.c 

fp.o: fp.c
	gcc -c -Wall -g fp.c


clean: 
	rm -f *.o
	echo Clean done