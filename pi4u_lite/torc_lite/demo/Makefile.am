bin_PROGRAMS= masterslave mbench1 fibo broadcast struct pipe async zerolength

AM_CFLAGS = @DEBUG_FLAG@ -I../include -I../uth -I.

LIBS = -L../src -ltorc $(MPILIB) -lpthread -lm 

masterslave_SOURCES = masterslave.c
masterslave_LDADD =$(LIBS)

mbench1_SOURCES = mbench1.c
mbench1_LDADD =$(LIBS)

fibo_SOURCES = fibo.c
fibo_LDADD =$(LIBS)

broadcast_SOURCES = broadcast.c
broadcast_LDADD =$(LIBS)

struct_SOURCES = struct.c
struct_LDADD =$(LIBS)

pipe_SOURCES = pipe.c
pipe_LDADD =$(LIBS)

async_SOURCES = async.c
async_LDADD =$(LIBS)

zerolength_SOURCES = zerolength.c
zerolength_LDADD =$(LIBS)

.c.o:
	$(CC) $(AM_CFLAGS) $(CFLAGS) -c $<
