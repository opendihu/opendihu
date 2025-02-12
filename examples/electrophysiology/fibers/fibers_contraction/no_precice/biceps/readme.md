# A biceps model

## Setup
- A realistic biceps geometry.
- Solver comprises both mechanics solver + fastmonodomain solver. 
- There are two options: one is using a shorten model and the other is using a hodgkin-huxley variant
- No preCICE involved. 

### How to build?
Follow OpenDiHu's documentation for installation, then run 
```
mkorn && sr
```
For a debug build, look into the documentation. 


### How to run?
To run the case go into the build directory and choose one of the two options:

- To run the shorten model:
```
./biceps_contraction ../settings_biceps_contraction.py ramp.py
```

- To run the hodgkin-huxley model:
```
./biceps_contraction_Fv ../settings_biceps_contraction.py ramp_Fv.py

```

If possible, this case should be run in parallel, as it is very time consuming. On the ipvs-epyc cluster a simulation with `end_time=50` ms using `mpirun -n 16 ./biceps_contraction ../settings_biceps_contraction.py ramp.py` ranks took around 1 hour, while using `mpirun -n 4 ./biceps_contraction ../settings_biceps_contraction.py ramp.py` took 4 hours approximately. The numer of ranks is limited by the mechanics mesh, which is very coarse. If you use 16 ranks, the ranks will split as n_x * n_y * n_z = 2 * 2 * 4. This splitting is not possible if you use the `left_biceps_brachii_7x7fibers.bin`, and you will have to either reduce the number of ranks or increase the mesh resolution. 


