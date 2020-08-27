gcc -g test_post.c -L../target/release -ls1_noisefloor -lpython3.6m -o test_post
LD_LIBRARY_PATH=../target/release valgrind --leak-check=full ./test_post

