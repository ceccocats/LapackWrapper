xavier test
===========
```
# install deps
sudo apt install libopenblas-dev 
# compile
make OPENBLAS=1 config
make
make tests

# run sample
./bin/test2-Timing_cuda
```

result
~~~~
Size N = 2
CPU MULT = 189.549 [ms] (lapack)
GPU MULT = 9053.26 [ms] (lapack)
All done!

Size N = 3
CPU MULT = 247.53 [ms] (lapack)
GPU MULT = 8867.59 [ms] (lapack)
All done!

Size N = 4
CPU MULT = 271.06 [ms] (lapack)
GPU MULT = 9030.32 [ms] (lapack)
All done!

Size N = 5
CPU MULT = 406.516 [ms] (lapack)
GPU MULT = 9056.74 [ms] (lapack)
All done!

Size N = 6
CPU MULT = 446.029 [ms] (lapack)
GPU MULT = 9122.81 [ms] (lapack)
All done!

Size N = 7
CPU MULT = 658.115 [ms] (lapack)
GPU MULT = 9221.58 [ms] (lapack)
All done!

Size N = 8
CPU MULT = 643.356 [ms] (lapack)
GPU MULT = 9272.12 [ms] (lapack)
All done!

Size N = 16
CPU MULT = 2766.27 [ms] (lapack)
GPU MULT = 9644.21 [ms] (lapack)
All done!

Size N = 32
CPU MULT = 16794.4 [ms] (lapack)
GPU MULT = 12100.1 [ms] (lapack)
All done!

Size N = 64
CPU MULT = 125060 [ms] (lapack)
GPU MULT = 20004 [ms] (lapack)
All done!


All done!
~~~~


submodule compilation
=====================

**LapackWrapper**

Prerequisite
------------

~~~~
cd third_parties
~~~~

**Linux**

~~~~
rake install_linux
~~~~

**OSX**

~~~
rake install_osx
~~~

**Windows**

* Visual Studio 2017 64 bit

~~~~
rake install_win[2017,x64]
~~~~

* Visual Studio 2017 32 bit

~~~~
rake install_win[2017,x86]
~~~~

COMPILE
=======

**On linux**

~~~~
make OPENBLAS=1 config
make clean
make
make install_local
~~~~

or using rake

~~~~
rake build_linux
~~~~

**On windows**

using MINGW on a bash shell

~~~~
make OPENBLAS=1 config
make clean
make
make install_local
~~~~

or using Visual Studio

~~~~
rake build_win[2017,x64]
rake build_win[2017,x86]
~~~~

**On OSX use**

~~~~
make ACCELERATE=1 config
make clean
make
make install_local
~~~~

**the library and header file are in the directory `lib`**

~~~~
lib/lib  (static)
   /dll  (dll windows import library)
         (unix shared library)
   /bin  (windows dynamic lib [the dll])
   /headers
~~~~
