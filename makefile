# Makefile 

F77        = gfortran
FFLAGS     = -c 
LDFLAGS    = 
SOURCES    = convert.f
OBJECTS    = $(SOURCES:.f=.o)
EXECUTABLE = convert

#all: $(SOURCE) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(F77) $(LDFLAGS) $(OBJECTS) -o $@

.f.o:
	$(F77) $(FFLAGS) $< -o $@

clean:
	rm -rf *o convert
