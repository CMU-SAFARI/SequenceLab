CC=g++
CFLAGS= -static-libstdc++ -static-libgcc -Wall -O3
#CFLAGS+= -fsanitize=address -g
LDFLAGS= -v -lz -lm -L libs/string_buffer/ -L src/ -lalign -lstrbuf -lpthread
EXEC=main
SRC=main.c src/edlib.cpp filters/base-counting/Base_Counting.c filters/SneakySnake/SneakySnake.c filters/qgram/qgram.c \
	filters/adjacency-filter/AdjacencyFilter.c filters/shd/SHD.c filters/magnet/MAGNET.c filters/shouji/Shouji.c \
	filters/hamming-distance/HD.c filters/magnet/MAGNET_DC.c filters/grim/grim.c filters/pigeonhole/pigeonhole.c \
	aligners/ksw2/ksw2_extd.c aligners/ksw2/ksw2_extd2_sse.c aligners/ksw2/ksw2_extz.c aligners/ksw2/ksw2_extz2_sse.c \
	aligners/ksw2/ksw2_gg.c aligners/ksw2/ksw2_gg2.c aligners/ksw2/ksw2_gg2_sse.c
OBJ=$(patsubst %.cpp, %.o, $(patsubst %.c, %.o, $(SRC)))

ifeq ($(sse2),)
	CFLAGS += -march=native
endif

ifneq ($(avx2),)
	CFLAGS += -mavx2
endif

all: $(EXEC)

main: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) $(CFLAGS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY: clean

clean:
	rm -f $(OBJ)

main.o: filters/base-counting/Base_Counting.h filters/SneakySnake/SneakySnake.h filters/qgram/qgram.h filters/adjacency-filter/AdjacencyFilter.h filters/shd/SHD.h filters/magnet/MAGNET.h src/edlib.h filters/shouji/Shouji.h filters/hamming-distance/HD.h filters/magnet/MAGNET_DC.h filters/grim/grim.h filters/pigeonhole/pigeonhole.h aligners/ksw2
kalloc.o: kalloc.h
ksw2_extd.o: ksw2.h
ksw2_extd2_sse.o: ksw2.h
ksw2_extf2_sse.o: ksw2.h
ksw2_extz.o: ksw2.h
ksw2_extz2_sse.o: ksw2.h
ksw2_gg.o: ksw2.h
ksw2_gg2.o: ksw2.h
ksw2_gg2_sse.o: ksw2.h
