# UpPipe

[![GitHub repository](https://img.shields.io/badge/GitHub-chi--0828%2FUpPipe-blue.svg)](https://github.com/chi-0828/UpPipe)
![GitHub top language](https://img.shields.io/github/languages/top/chi-0828/UpPipe?color=blue&logo=Ionic&logoColor=white)
![GitHub commit activity (branch)](https://img.shields.io/github/commit-activity/w/chi-0828/UpPipe)
![GitHub last commit (by committer)](https://img.shields.io/github/last-commit/chi-0828/UpPipe)
[![C++ version](https://img.shields.io/badge/c++-14-yellow)](https://docs.npmjs.com/)
[![g++ version](https://img.shields.io/badge/gcc-8.3.0-yellow)](https://docs.npmjs.com/)
<br>
UpPipe is an RNA abundance quantification design on a real processing-near-memory system ([UPMEM DPU](https://www.upmem.com/)); the paper of this project is published in [Design Automation Conference (DAC) 2023](https://www.dac.com/)

## Citation
> Liang-Chi Chen,  Chien-Chung Ho, and Yuan-Hao Chang, “UpPipe: A Novel Pipeline Management on In-Memory Processors for RNA-seq Quantification," ACM/IEEE Design Automation Conference (DAC), San Francisco, CA, USA, July 9-13, 2023.
- Online proceedings of DAC'23 are available soon

## Materials
- Link to paper (coming soon)
- [Link to slides](https://drive.google.com/file/d/1XaUErirVkLod5UZwsReGUwLDN2Af026Q/view?usp=drive_link)
- [Link to poster](https://drive.google.com/file/d/1OGtMobOE1xZWm_qes1gTFDT9nAnk1r31/view?usp=drive_link)
- [DPU-based kallisto](https://github.com/chi-0828/RNA-Abundance-Quantification-on-UPMEM)
- [UPMEM use cases](https://www.upmem.com/ressources/)
- [UPMEM SDK](https://sdk.upmem.com/2021.3.0/index.html)

## Hardware/System Prerequisites
The project has to be run on a system equipped with UPMEM DRAM Processing Units (DPUs), and the kernel system requires installing the [UPMEM SDK](https://sdk.upmem.com/2021.3.0/index.html)

## Start
```=shell
git clone https://github.com/chi-0828/UpPipe.git
cd UpPipe
chmod +x build.sh
./build.sh
make -j4
```

## Usage
### Allocate transcriptome to DPU(s)
- `KMER SIZE` should be 3, 5, ..., 31
- `NUMBER OF DPU(s) in a PIPELINE WORKER` should be less than 64 in our suggestion
```=shell
./UpPipe build \
            -k KMER SIZE  \
            -i OUTPUT INDEX FILE PATH \
            -d NUMBER OF DPU(s) in a PIPELINE WORKER \
            -f TRANSCRIPTOME FILE PATH
```
### Run alignment step for quantification
- The size of k-mer is already set in `INPUT INDEX FILE`, this setting cannot be changed in this step 
```=shell
./UpPipe alignment \
            -i INPUT INDEX FILE PATH \
            -r NUMBER OF PIPELINE WORKER(s) \
            -f INPUT RNA READ FILE PATH
```
### Parameters setting (dpu_app/dpu_def.h)
- `KMER SIZE` less than 7 may lead to inaccurate mapping result
- `NUMBER OF DPU(s) in a PIPELINE WORKER` should be less than 64 for optimal performance
- The `number of transcript / NUMBER OF DPU(s) in a PIPELINE WORKER` must be less than 200 (`COUNT_LEN in dpu_app/dpu_def.h`)
- Setting `READ_LEN` to the sequence length of RNA reads
- Setting `WRAM_READ_LEN` to the a number which is larger than `READ_LEN` and divisible by 8
- `WRAM_PREFETCH_SIZE` is the size for WRAM pre-feteching, 16 is the optimal size in most situations

## Test
- To build the index file by 11-mer and allocate to 60 DPUs
```=shell
./UpPipe build \
            -k 11  \
            -i test/test.idx \
            -d 60 \
            -f test/tran.fa
```
- To run alignment with 40 pipeline workers
```=shell
./UpPipe alignment \
            -i test/test.idx \
            -r 40 \
            -f test/read.fa
```
- Performance: UpPiep uses 40 pipeline workers
```
real    0m2.747s
```
- Performance: UpPiep uses 20 pipeline workers
```
real    0m3.584s
```
- Performance: [kallisto](https://github.com/pachterlab/kallisto)
```
real    0m4.003s
```
- To note that UpPipe shows its efficiency more in the large size dataset due to the porcessing-in-memory features

