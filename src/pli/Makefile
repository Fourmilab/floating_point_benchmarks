
BASE    = $(HOME)/linuxtools/pli/pli-0.9.9b
PLI     = $(BASE)/plic
PLIFLGS = -lsiaxgo -ew "-cn(^) -i$(BASE)/lib/include"
INC     = $(BASE)/lib/include
LIBDIR  = $(BASE)/lib/static

PROGRAMS = fbench

all:    $(PROGRAMS)

%.o:    %.pli
	${PLI} -C -dELF ${PLIFLGS}  $^ -o $*.o

%:      %.o
	ld -z muldefs -Bstatic -M -o $@         \
	--oformat=elf32-i386                    \
	-melf_i386                              \
	-e main                                 \
	 $@.o                                   \
	--start-group                           \
	  $(BASE)/lib/libprf.a                  \
	--end-group                             \
	>$@.map

clean:
	rm -f $(PROGRAMS) *.lst *.map

time:   $(PROGRAMS)
	time -p ./fbench
