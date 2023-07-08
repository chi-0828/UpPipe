INCLUDE := -I.
LDFLAGS := 
CFLAGS := -O3
OBJ := ../../obj/dpu_app/
target = $(OBJ)dpu
CC := dpu-upmem-dpurte-clang -DNR_TASKLETS=$(TASKLETS)

$(target): $(OBJ)hash.o $(OBJ)dpu.o main.c
	$(CC) -o $(target) $^ $(LDFLAGS) $(CFLAGS)
$(OBJ)hash.o : hash.c
	$(CC) -o $(OBJ)hash.o -c hash.c $(CFLAGS) $(INCLUDE)
$(OBJ)dpu.o : dpu.c
	$(CC) -o $(OBJ)dpu.o -c dpu.c $(CFLAGS) $(INCLUDE)