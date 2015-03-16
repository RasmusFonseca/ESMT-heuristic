CC=g++
# Add -DNDEBUG and -O2 in CFLAGS before testing
CFLAGS=-Wall -I./ -Wfatal-errors -O2 -DNDEBUG

MAIN=main

OBJECTS = $(subst cpp,o,$(wildcard *.cpp))
OBJECTS+= $(subst cpp,o,$(wildcard steiner/*.cpp))
OBJECTS+= $(subst cpp,o,$(wildcard steiner/*/*.cpp))
OBJECTS+= $(subst cpp,o,$(wildcard test/*.cpp))

all: ${MAIN} $(OBJECTS)

%.o: %.cpp
	g++ $(CFLAGS) -c -o $@ $^

$(MAIN): $(OBJECTS)
	$(CC) -o main $(OBJECTS)

clean:
	rm -f $(MAIN) $(OBJECTS)

run:
	./$(MAIN)
