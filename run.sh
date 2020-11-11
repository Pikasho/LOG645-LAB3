make all
echo ""
echo "<<<"
mpirun -np 9 ./lab3 9 9 3 0.00025 0.1
echo ">>>"

echo ""
echo "<<<"
mpirun -np 15 ./lab3 15 9 10 0.00025 1
echo ">>>"
make clean

# make all
# echo "mpirun -np 48 ./lab3 380 420 6 0.00025 0.1"
# echo "<<<"
# mpirun -np 48 ./lab3 380 420 6 0.00025 0.1 >> test1_v1.txt
# echo ">>>"

# echo "mpirun -np 6 ./lab3 25 19 700 0.00025 0.1"
# echo "<<<"
# mpirun -np 6 ./lab3 25 19 700 0.00025 0.1 >> test2_v1.txt
# echo ">>>"

# echo "mpirun -np 41 ./lab3 1000 15 15 0.00025 0.1"
# echo "<<<"
# mpirun -np 64 ./lab3 1000 15 15 0.00025 0.1 >> test3_v1.txt 
# echo ">>>"

# echo "mpirun -np 35 ./lab3 10 1200 35 0.00025 0.1"
# echo "<<<"
# mpirun -np 35 ./lab3 10 1200 35 0.00025 0.1 >> test4_v1.txt
# echo ">>>"
# make clean