CFLAGS=-Wall -g
LDFLAGS=-lR -lm
MEMFLAGS=-fsanitize=address -static-libasan

INCPATH=-I/usr/share/R/include
INCDIR=-Iinclude
OBJDIR=obj
SRCDIR=src
BINDIR=bin

CC=gcc -fopenmp
EXEC=clean dir $(BINDIR)/test
SRC:=$(wildcard $(SRCDIR)/*.c)
OBJ:=$(SRC:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

all: $(EXEC)

#run: $(BINDIR)/test
#	./$(BINDIR)/test ../Exemples/ADENOS10.xyz 3 13

$(BINDIR)/test: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

$(OBJ) : $(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(INCPATH) $(INCDIR) -o $@ -c $< $(CFLAGS)

dir:
	mkdir $(OBJDIR)
	mkdir $(BINDIR)

.PHONY: clean mrproper all mem

clean:
	rm -rf $(OBJDIR)
	rm -rf $(BINDIR)
	ls -s

mrproper: clean
	rm -rf results