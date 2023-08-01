target := UpPipe
SRC_DIR := src
OBJ_DIR := obj
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
INCLUDE := -I./ext/htslib/ -I.
LDFLAGS := `dpu-pkg-config --cflags --libs dpu` ./ext/htslib/libhts.a -pthread -lz -lrt -g -O3
CXXFLAGS := -g `dpu-pkg-config --cflags --libs dpu` -O3

DPU := src/dpu_app/dpu

all : $(target) $(DPU)

$(target): $(OBJ_FILES)
	g++ -o $(target) $^ $(LDFLAGS) 

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp 
	g++ -c -o $@ $< $(CXXFLAGS) $(INCLUDE)

$(DPU):
	$(MAKE) -C src/dpu_app/ TASKLETS=1

clean: 
	rm -f $(target) obj/*.o obj/dpu_app/*