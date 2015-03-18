CC=g++
# Add -DNDEBUG and -O2 in CFLAGS before testing
CFLAGS=-Wall -I./ -Wfatal-errors -O2 -DNDEBUG

MAIN=esmt-heuristic

OBJECTS = $(subst cpp,o,$(wildcard *.cpp))
OBJECTS+= $(subst cpp,o,$(wildcard steiner/*.cpp))
OBJECTS+= $(subst cpp,o,$(wildcard steiner/*/*.cpp))
OBJECTS+= $(subst cpp,o,$(wildcard test/*.cpp))

CFLAGS+= -std=c++11


all: $(MAIN) $(OBJECTS)

%.o: %.cpp
	g++ $(CFLAGS) -c -o $@ $^

$(MAIN): $(OBJECTS)
	$(CC) -o $(MAIN) $(OBJECTS)

clean:
	rm -f $(MAIN) $(OBJECTS)

run:
	./$(MAIN)
