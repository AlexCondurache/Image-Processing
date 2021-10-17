main_iplib: bmp.o ip_lib.o main_iplib.o
	gcc bmp.o ip_lib.o main_iplib.o -o main_iplib -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

bmp.o: bmp.c
	gcc bmp.c -o bmp.o -c -Wall
	
main_iplib.o: main_iplib.c
	gcc main_iplib.c -o main_iplib.o -c -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
	
ip_lib.o: ip_lib.c
	gcc ip_lib.c -o ip_lib.o -c -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
	
clean:
	@rm -f main_iplib bmp.o ip_lib.o main_iplib.o
