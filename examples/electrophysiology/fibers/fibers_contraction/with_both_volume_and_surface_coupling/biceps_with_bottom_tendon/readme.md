# Biceps with bottom tendon

## Setup
We use preCICE to couple three participants (muscle mechanics, muscle fibers, bottom tendon). The muscle and the tendon use the realistic geometries obtained from imaging and post-processed by Benjamin Maier.  

The `precice_config.xml` uses a composition of coupling scheme: we combine a `serial-explicit` scheme with a `parallel-implicit` scheme.

### How to build?
Follow OpenDiHu's documentation for installation, then run 
```
mkorn && sr
```

### How to run?
You will need one terminal per each participant. All terminals should be at the same directory, e.g., `biceps_with_bottom_tendon/build_release`.

terminal 1: muscle fibers
```
./partitioned_fibers ../settings_partitioned_fibers.py ramp.py
```
terminal 2: muscle mechanics
```
./muscle_contraction ../settings_muscle_contraction.py ramp.py
```
terminal 3: bottom tendon
```
./tendon ../settings_tendon_implicit_dirichlet_neumann.py tendon.py
```

### How to run in a cluster?
Here's an example on how to lunch the three participants in the cluster (Add this code to your bash file):
```
echo "Launching muscle fibers"
mpirun -n 16 ./partitioned_fibers ../settings_partitioned_fibers.py ramp.py &> muscle_fibers.log &

echo "Launching muscle mechanics"
mpirun -n 1 ./muscle_contraction ../settings_muscle_contraction.py ramp.py &> muscle_mechanics.log &

echo "Launching bottom tendon"
mpirun -n 1 ./tendon ../settings_tendon_implicit_dirichlet_neumann.py tendon.py
 &> tendon.log

echo "Simulation completed."

```
