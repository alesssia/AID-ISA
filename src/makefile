OBJS = Matrix.o Cluster.o Bicluster.o Driver.o AID-ISA.o

CFLAGS = -g -Wall -O2 -Wno-unused-function
LIBS = -lboost_program_options 

all: aid_isa

aid_isa: $(OBJS)
	/bin/rm -rf ../bin/
	mkdir ../bin
	g++ $(CFLAGS) $(LIBS) -o AID-ISA $(OBJS)
	mv AID-ISA ../bin/

%.o: %.cpp
	g++ $(CFLAGS) -c $<

clean:
	/bin/rm -f $(OBJS) utilities.h.gch
	/bin/rm -rf ../bin/
