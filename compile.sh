icc -o in2mart in2mart.c -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread  -mcmodel=medium
icc -o i2mart i2mart.c -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread  -mcmodel=medium
icc -o i2jmart i2jmart.c -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread  -mcmodel=medium
rm -f *.c~
rm -f *.t3c~
