#! /bin/sh
./vseg.sh

make clean

make

chmod +x ./test1.sh
chmod +x ./test2.sh
chmod +x ./test3.sh
chmod +x ./test4.sh

./test1.sh &
./test2.sh &
./test3.sh &
./test4.sh &