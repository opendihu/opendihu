# Biceps with tendons

## Setup
We use preCICE to couple 5 participants (muscle mechanics, muscle fibers, bottom tendon, top-a tendon, top-b tendon). The muscle and the tendons use the realistic geometries obtained from imaging and post-processed by Benjamin Maier.  

This case is equivalent to the [multiple tendons with electrophysiology case](https://github.com/opendihu/opendihu/tree/develop/examples/electrophysiology/fibers/fibers_contraction/with_tendons_precice/multiple_tendons_with_electrophysiology), but with the addition of volume coupling within the muscle participant. 

The `precice_config.xml` uses a composition of coupling scheme: we combine a `serial-explicit` scheme with a `multi` scheme.

### How to build?
Follow OpenDiHu's documentation for installation, then run 
```
mkorn && sr
```

### How to run?
You will need one terminal per each participant. All terminals should be at the same directory, e.g., `biceps_with_bottom_tendon/build_release`.

terminal 1: muscle fibers
```
./muscle_fibers ../settings_muscle_fibers.py ramp.py
```
terminal 2: muscle mechanics
```
./muscle_mechanics ../settings_muscle_mechanics.py ramp.py
```
terminal 3: bottom tendon
```
./tendon_linear ../settings_tendon_bottom.py ramp.py
```
terminal 4: top-a tendon
```
./tendon_linear ../settings_tendon_top_a.py ramp.py
```
terminal 5: top-b tendon
```
./tendon_linear ../settings_tendon_top_b.py ramp.py
```

### How to run in a cluster?
Here's an example on how to lunch the three participants in the cluster (Add this code to your bash file):
```
echo "Launching left muscle"
mpirun -n 16 ./partitioned_fibers ../settings_partitioned_fibers.py ramp.py &> muscle_fibers.log &

echo "Launching tendon"
mpirun -n 1 ./muscle_contraction ../settings_muscle_contraction.py ramp.py &> muscle_mechanics.log &

echo "Launching right muscle"
mpirun -n 1 ./tendon ../settings_tendon_implicit_dirichlet_neumann.py tendon.py
 &> tendon.log

echo "Simulation completed."

```
