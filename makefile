# project name (generate executable with this name)
TARGET   = cageMol.exe

CFLAGS=-Wall -g
LDFLAGS=-lR -lm

INCPATH=-I/usr/share/R/include
INCDIR=-Iinclude
OBJDIR=obj
SRCDIR=src
BINDIR=bin

CC=gcc -fopenmp -O3
EXEC=clean dir $(BINDIR)/$(TARGET)
SRC:=$(wildcard $(SRCDIR)/*.c)
OBJ:=$(SRC:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

all: $(EXEC)

demo: all $(BINDIR)/$(TARGET)
	./$(BINDIR)/$(TARGET) -i ./demos/substrates/ADENOS10.xyz

$(BINDIR)/$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

$(OBJ) : $(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(INCPATH) $(INCDIR) -o $@ -c $< $(CFLAGS)

dir:
	mkdir $(OBJDIR)
	mkdir $(BINDIR)

.PHONY: clean mrproper all

pathfinding:  #-DDEBUGAstar -DDEBUGDij
	$(CC) $(INCPATH) $(INCDIR) -o test/initialization.o -c src/initialization.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/input.o -c src/input.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/output.o -c src/output.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/pathFinding.o -c src/pathFinding.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/structure.o -c src/structure.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/structureAsp.o -c src/structureAsp.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/structureGph.o -c src/structureGph.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/structureLst.o -c src/structureLst.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/structureMN.o -c src/structureMN.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/structureMol.o -c src/structureMol.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/structureNH.o -c src/structureNH.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/structurePT.o -c src/structurePT.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/structureShl.o -c src/structureShl.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/util.o -c src/util.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/voxelization.o -c src/voxelization.c $(CFLAGS)
	$(CC) $(INCPATH) $(INCDIR) -o test/pathFindingTest.o -c test/pathFindingTest.c $(CFLAGS)
	$(CC) -o test/pathFinding.exe test/initialization.o test/input.o test/output.o test/pathFinding.o test/structure.o test/structureAsp.o test/structureGph.o test/structureLst.o test/structureMN.o test/structureMol.o test/structureNH.o test/structurePT.o test/structureShl.o test/util.o test/voxelization.o test/pathFindingTest.o $(LDFLAGS)
	./test/pathFinding.exe
	rm test/*.o

clean:
	rm -rf $(OBJDIR)
	rm -rf $(BINDIR)
	ls

mrproper: clean
	rm -rf results