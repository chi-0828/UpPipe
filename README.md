# UpPipe: In-Memory Processing based RNA-seq Quantification
[![GitHub repository](https://img.shields.io/badge/GitHub-chi--0828%2FUpPipe-blue.svg)](https://github.com/chi-0828/UpPipe)
![GitHub top language](https://img.shields.io/github/languages/top/chi-0828/UpPipe?color=blue&logo=Ionic&logoColor=white)
![GitHub commit activity (branch)](https://img.shields.io/github/commit-activity/w/chi-0828/UpPipe)
![GitHub last commit (by committer)](https://img.shields.io/github/last-commit/chi-0828/UpPipe)
[![C++ version](https://img.shields.io/badge/c++-14-yellow)](https://docs.npmjs.com/)
[![g++ version](https://img.shields.io/badge/gcc-8.3.0-yellow)](https://docs.npmjs.com/)
<br>
UpPipe is a DPU-based RNA abundance quantification design; the paper of this project is published in [Design Automation Conference (DAC) 2023](https://www.dac.com/)

## Citation
- Online proceedings of DAC'23 are available soon

## Materials
- [Link to paper]()
- [Link to slides](https://drive.google.com/file/d/1XaUErirVkLod5UZwsReGUwLDN2Af026Q/view?usp=drive_link)
- [Link to poster](https://drive.google.com/file/d/1OGtMobOE1xZWm_qes1gTFDT9nAnk1r31/view?usp=drive_link)
- [UPMEM SDK](https://sdk.upmem.com/2021.3.0/index.html)
- [DPU-based kallisto](https://github.com/chi-0828/RNA-Abundance-Quantification-on-UPMEM)

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
### Suggestion
- `KMER SIZE` less than 7 may lead to inaccurate mapping result
- `NUMBER OF DPU(s) in a PIPELINE WORKER` should be less than 64 in optimal situations
- However, the `number of transcript / NUMBER OF DPU(s) in a PIPELINE WORKER` must be less than 104 (`T_LEN in dpu_app/dpu_def.h`), so `NUMBER OF DPU(s) in a PIPELINE WORKER` may exceed 64 in some cases
- Making `T_LEN` greater than 104 may cause the DPU to fault due to insufficient WRAM
- The project is still in progress for extension work

## Test
- To build the index file by 7-mer and allocate to 60 DPUs
```=shell
./UpPipe build \
            -k 13  \
            -i test/test.idx \
            -d 60 \
            -f test/tran.fa
```
- To run alignment with 10 pipeline workers
```=shell
./UpPipe alignment \
            -i test/test.idx \
            -r 10 \
            -f test/read.fa
```

