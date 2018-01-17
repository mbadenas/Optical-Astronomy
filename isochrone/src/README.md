# pyIsochrone
### Compilation instructions with f2py

To compile the module **pyIsochrone** (`.so`) file, execute the following line inside the folder `isochrone/src/`

    f2py -c pyiso.pyf pymodels.f

Afterwards, move the resulting `pyIsochrone.so` to the `isochrone/` folder in order to allow the subroutine **pymodels**
to access the files in `isochrone/data/`.

Alternatively, this same f2py command can be written directly in the folder `isochrone/` or any other folder specifying the 
path to `pyiso.pyf` and `pymodels.f` when running the command.
