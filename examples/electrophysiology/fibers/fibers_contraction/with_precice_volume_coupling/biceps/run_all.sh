cd build_release
rm -r precice-run
echo "Launching mechamics"
mpirun -n 1 ./muscle_contraction ../settings_muscle_contraction.py ramp.py &>mechanics.log &
echo "Launching fibers"
mpirun -n 1 ./partitioned_fibers ../settings_partitioned_fibers.py ramp.py
echo "Simulation completed."

