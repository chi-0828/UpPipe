mkdir obj

cd ./ext/htslib
autoheader    
autoconf       
./configure    
make
make install
cd ..