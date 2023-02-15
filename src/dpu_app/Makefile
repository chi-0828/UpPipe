target = dpu
INCLUDE := -I.
LDFLAGS := 
CFLAGS := -O3
CC := dpu-upmem-dpurte-clang -DNR_TASKLETS=$(TASKLETS)

$(target): hash.o dpu.o main.c
	$(CC) -o $(target) $^ $(LDFLAGS) $(CFLAGS)
hash.o : hash.c
	$(CC) -c hash.c $(CFLAGS) $(INCLUDE)
dpu.o : dpu.c
	$(CC) -c dpu.c $(CFLAGS) $(INCLUDE)

clean: 
	rm -f dpu *.o