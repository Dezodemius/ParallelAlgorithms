mkdir bin

echo "COMPILING"

if ! gcc main.c -o bin/main.out
then
  exit
fi

echo "DONE"

echo "WITHOUT OPTIMIZATION RUN"
./bin/main.out -O0
echo "DONE"

echo "-O1 OPTIMIZATION RUN"
./bin/main.out -O1
echo "DONE"

echo "-O2 OPTIMIZATION RUN"
./bin/main.out -O2
echo "DONE"

echo "-O3 OPTIMIZATION RUN"
./bin/main.out -O3
echo "DONE"
