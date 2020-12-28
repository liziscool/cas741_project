SOURCE_MAIN = transform.c plot.c main.c 
#INCLUDES = $(SOURCES_MAIN:.c=.h)
INCLUDES = ../include/transform.h  ../include/plot.h
CC = gcc
#CXX = g++
CFLAGS  = -fPIC  -Wall -O3  -g 
LINK = gcc
LDFLAGS = 
LIB =  -lfftw3 

INCLUDE_DIR =  -I.

OBJECT_PATH = objects
vpath .. %.o $(OBJECT_PATH)

OBJECTS_MAIN = $(patsubst %.c,objects/%.o,$(SOURCE_MAIN)) 

TARGET =   run

$(TARGET): $(OBJECTS_MAIN)
	$(LINK) $(LDFLAGS)  -o $(TARGET) $(LIB)  $(OBJECTS_MAIN) $(PNG)


all: $(TARGET) 

.SUFFIXES: .h .c

objects/%.o:  %.c
	$(CC)  $(INCLUDE_DIR)  $(CFLAGS) -c $(<) -o $@ -DNO_FREETYPE

clean:
	rm -f objects/*.o $(TARGET)

backup:
	ci -l *.c *.h *.inc show.html makefile

depend:	
	makedepend -I. -Y.  -pobjects/ -fmakefile $(INCLUDES) $(SOURCE_MAIN)

clean_png:
	find . -type f -name '*.png' -delete


# DO NOT DELETE
objects/main.o: ../include/transform.h  ../include/plot.h
objects/transform.o: ../include/plot.h 
objects/plot.o:




