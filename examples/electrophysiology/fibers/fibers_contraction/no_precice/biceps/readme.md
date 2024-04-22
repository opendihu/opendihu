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

