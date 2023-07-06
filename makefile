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

clean:
	rm -rf $(OBJDIR)
	rm -rf $(BINDIR)
	ls

mrproper: clean
	rm -rf results