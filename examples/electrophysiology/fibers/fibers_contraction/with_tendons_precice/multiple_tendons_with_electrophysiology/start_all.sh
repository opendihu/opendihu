rm -r precice-run
echo "Launching muscle"
mpirun -n 1 ./muscle_electrophysiology_precice settings_muscle.py ramp.py --case_name "default"  &> muscle.log &
echo "Tendon 1"
mpirun -n 1 ./tendon_linear_precice_dynamic settings_tendon_bottom.py --case_name "default"  &> tendon_bottom.log &
echo "Tendon 2a"
mpirun -n 1 ./tendon_linear_precice_dynamic settings_tendon_top_a.py --case_name "default" &> tendon_top_a.log &
echo "Tendon 2b"
mpirun -n 1 ./tendon_linear_precice_dynamic settings_tendon_top_b.py --case_name "default" &> tendon_top_b.log
echo "Simulation completed."
