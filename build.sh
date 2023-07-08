mkdir obj
mkdir obj/dpu_app
cd ./ext/htslib
autoheader    
autoconf       
./configure    
make
make install
cd ..