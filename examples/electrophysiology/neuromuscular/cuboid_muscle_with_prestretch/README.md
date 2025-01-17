This folder contains a cuboid muscle example. We apply a prestretch and a contraction solver. The prestretch and the contraction solver have both the same structure: They couple the `FastMonodomainSolver`+ `MuscleContractionSolver`. In the prestretch one we set `dynamic = False`, and in the contration one we set `dynamic = True`. The boundary conditions are also different for the prestretch and the contraction solver. In the prestretch part of the simulation, we fix one side and pull the other with a constant force. In the contraction part of the simulation, we remove the pulling force.

To compile the code:

```
mkorn && sr
```

To run the code: 
```
cd build_release
./muscle_with_prestretch ../settings_muscle_with_prestretch.py --force 1.1
```

'--force' specifies the pulling force. A useful range of values would be from 1 to 10 N. 
A requirement to run this case is to have set up `OPENDIHU_HOME` in your system. 