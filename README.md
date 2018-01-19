# SIKE_generalized

Supersingular isogeny graph key exchange, using prime numbers (2,3) as in the spec, or any odd prime, as in the paper from Costello and Hisil.

### Run instruction

Contains a makefile.
To run the protocole, type
$ make
$./Main

To run the protocole with the profiler, use "makefile_with_profiling" instead.

### Content

Header.h with declarations.

fp.c to initialize, handle and clear elements in fp2.

curve_point.c to initialize, handle and clear curve and points.

montgomery.c contains methods for Montgomery's arithmetic.

isogeny.c contains methods to compute images from curves and points by 2, 3? or any odd degree isogeny.

pk_sk_param.c  contains methods to initialize, set and clear public parameters, private key and public key, including key generation and key exchange.

Main.c contains wrappers for key generation and key exchange, and a set of tests to run the protocole.

Makefile 
