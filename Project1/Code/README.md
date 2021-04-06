# FYS4411 VMC Code

This VMC code is built on the code skeleton made by Morten Ledum: https://github.com/mortele/variational-monte-carlo-fys4411.

The typical flow of the code is as follows:

You run main(`./vmc`), where a VMC system is assembled. A `System` object is created, which then takes a `Hamiltionan`, `InitialState`, and `WaveFunction` object, and some more parameters for the system. Then, the `System` is handed to a `ParamTester` object, which "runs" the `System` for the correct values of the variational parameter alpha. By "running" the system, we mean running the metropolis algorithm, which is handeled entirely within the `System` class itself. From there it calls on functions from the objects it got in main to calculate the values it needs to perform the metropolis algorithm and sample various expected values. The sampling is handled by a `Sampler` object, which also returns the gradient to the `ParamTester`.

We have parallelized the metropolis algorithm using OpenMP. Each thread has its own particle objects(whose position is what defines the state of the quantum system), but they share most other values and objects. When a thread needs to store a value, we use vectors or arrays in shared objects like the `Sampler` or correlated wave function where each thread has its own index.

Every time `ParamTester` finishes a set of "runs" of the system, an output file is created in the folder "Output". `alphaGrid()` creates a file with the ending "_Grid" (a grid of alpha values and energies), `alphaGD()` makes the ending "_GD" (the alphas during gradient descent), and `bigCalc()` makes a file with the ending "_Big"(all energy samples from a large calculation) and "_oB"(the onebody density samples, which make a histogram). These files are read in the various Jupyter Notebooks in the folder "Notebooks". Utility functions from the file "customRead.py" are used to extract results from the files easily.

Currently, when you run main, only a small demo calculation is run. To reproduce our results, you must call one of the many other functions defined in main.

### Compiling the project using CMake
Run the script compile_project via `./compile_project`. Now the project can be run by executing `./vmc` in the top-directory.

### Unit tests
We have written unit tests in the folder UnitTests. The code there can be compiled, run and cleaned in the same way as in the top-directory. Simply run the script compile_project via `./compile_project` in the folder UnitTests. Then you can run the tests by executing `./unitTests`.

#### Cleaning the directory
Run `make clean` in the top-directory to remove the executable `vmc` and the `build`-directory.

#### Windows
Compilation of the project using Windows is still an open question to us. CMake should be OS-independent, but `make` does not work on Windows.
