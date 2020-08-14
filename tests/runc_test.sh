gcc test_clib.c -L../target/release -ls1_noisefloor -lpython3.6m -o test_clib
LD_LIBRARY_PATH=../target/release valgrind ./test_clib
