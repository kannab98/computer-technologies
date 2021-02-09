CC=g++
CFLAGS=
LDFLAGS=
SOURCES=newton.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=$(SOURCES:.cpp=)

all: 
	$(CC) $(SOURCES) -o $(EXECUTABLE)

#$(EXECUTABLE): $(OBJECTS)
	#$(CC) $(LDFLAGS) $(OBJECTS) -o $@

#.cpp.o:
	#$(CC) $(CFLAGS) $< -o $@

	
