# wscad2016-short-course

advisor

```
advixe-gui
```

Memory Access Optimization

Different strides on the same code (advisor: coordinate-aos-soa)

```
icc stride-test.c -o stride-test -g
```

optimization with O3

```
icc autoO3.c -o autoO3 -O0

time ./autoO3

real    0m8.428s
user    0m8.422s
sys     0m0.001s
```

optimization with O3 and vec-report 

```
icc autoO3.c -o autoO3 -O3 -vec-report6

time ./autoO3

real    0m0.003s
user    0m0.001s
sys     0m0.001s
```

```
cat autoO3.optrpt

Begin optimization report for: main()

    Report from: Vector optimizations [vec]


LOOP BEGIN at autoO3.c(10,3)
   remark #15344: loop was not vectorized: vector dependence prevents vectorization
   remark #15346: vector dependence: assumed OUTPUT dependence between A line 14 and A line 14
   remark #15346: vector dependence: assumed OUTPUT dependence between A line 14 and A line 14

   LOOP BEGIN at autoO3.c(13,5)
      remark #15414: loop was not vectorized: nothing to vectorize since loop body became empty after optimizations
   LOOP END
LOOP END
```

auto-vectorization (advisor: auto-vectorization)

```
icc autovec.c -o autovec -vec-report=5
cat autovec.optrpt

icc autovec.c -o autovec -vec-report=5 -O3
cat autovec.optrpt

icc autovec.c -o autovec -vec-report=5 -g -O3
cat autovec.optrpt

icc autovec.c -o autovec -vec-report=5 -g -O3 -xhost
cat autovec.optrpt

gcc autovec.c -O3 -o autovecGCC  
gcc autovec.c -O3 -ftree-vectorizer-verbose=1 -o autovecGCC 
```

Examples of loops that do not auto vectorize

```
icc novec.c -o novec -vec-report=6 -g -O3
cat novec.optrpt
```

memory aliasing
```
icc func.c -c func.o -vec-report=6 -g -O3
cat func.optrpt
icc func.c -c func.o -vec-report=6 -g -O3 -fargument-noalias
cat func.optrpt
include restrict on quad
icc func.c -c func.o -vec-report=6 -g -O3 -restrict
cat func.optrpt
```

ivdep
```
include ivdep on addfloats
icc func.c -c func.o -vec-report=6 -g -O3 -restrict
cat func.optrpt
```

Guided Vectorization

Matrix

```
make icc
check dependencies
pragma simd on inner loop

```

Interpolation 
```
  make clean ; make OMP=no
  cat main.optrpt 
  time ./main 
  
  make clean ; make OMP=no
  cat main.optrpt 
  time ./main 
```

Particle Binning (advisor: particle-binning)

```
change to reference version
make clean ; make app-CPU
time ./app-CPU 

change to vectorized version
make clean ; make app-CPU
time ./app-CPU 
```

Diffusion

Lattice Boltzmann


References

https://software.intel.com/en-us/articles/memory-management-for-optimal-performance-on-intel-xeon-phi-coprocessor-alignment-and

https://software.intel.com/en-us/articles/data-alignment-to-assist-vectorization

https://software.intel.com/en-us/articles/memory-management-for-optimal-performance-on-intel-xeon-phi-coprocessor-alignment-and

http://www.cs.uu.nl/docs/vakken/mov/
