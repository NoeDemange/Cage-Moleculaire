CFLAGS=-Wall -g
LDFLAGS=-lR -lm
MEMFLAGS=-fsanitize=address -static-libasan

INCPATH=-I/usr/share/R/include
INCDIR=-I../include
OBJDIR=../obj
SRCDIR=.

CC=gcc -fopenmp
EXEC=clean dir test
SRC:=$(wildcard $(SRCDIR)/*.c)
OBJ:=$(SRC:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

all: $(EXEC)

test: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

$(OBJ) : $(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(INCPATH) $(INCDIR) -o $@ -c $< $(CFLAGS)

dir:
	mkdir $(OBJDIR)

.PHONY: clean mrproper all mem

clean:
	rm -rf $(OBJDIR)
	ls -s

mrproper: clean
	rm -rf $(EXEC) ../Results