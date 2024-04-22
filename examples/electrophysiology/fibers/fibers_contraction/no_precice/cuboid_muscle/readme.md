# A cuboid muscle model

## Setup
- A dummy cuboid muscle geometry. 
- Solver comprises both mechanics solver + fastmonodomain solver. 
- It uses the CellML model "hodgkin_huxley-razumova".
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
./muscle_contraction ../settings_muscle_contraction.py variables.py
```

> [!WARNING]  
> Currently fails to run in parallel. 
