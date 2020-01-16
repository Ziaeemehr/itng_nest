[link](https://github.com/nest/nestml)
[extension_modules](https://nest.github.io/nest-simulator/extension_modules)

### Adding paths to the .bashrc file
This is the most important step. With out adding the path you get lots of 
error saying the header files did not find.
suppose you have installed nest here:
-  `/home/username/path_to_your/nest-simulator-2.18.0_build`

-  in `.bashrc` or `.zshrc`:

```
export NEST_INSTALL_DIR=/home/username/path_to_your/nest-simulator-2.18.0_build
export NEST_MODULE_PATH=//home/username/path_to_your/nest-simulator-2.18.0_build/lib/nest:$NEST_MODULE_PATH
export SLI_PATH=/home/username/path_to_your/nest-simulator-2.18.0_build/share/nest/sli:$SLI_PATH
export CPLUS_INCLUDE_PATH=//home/username/path_to_your/nest-simulator-2.18.0_build/include/nest
```

### Prerequisites
Download, build and install NEST. NEST should be built outside the source code directory.
Install CMake version 2.8.12 or later.
The NEST source code and installation directory must be accessible for building modules.

### Building MyModule
As a starting point, try to build MyModule as follows:

1. From the NEST source directory, copy directory examples/MyModule to somewhere outside the NEST source, build or install directories.

2. Create a build directory for it on the same level as MyModule (e.g. mmb)

```{r, engine='bash', count_lines}
cd /path/to/MyModule
cd ..
mkdir mmb
cd mmb
```


3. Configure. The configure process uses the script nest-config to find out where NEST is installed, where the source code resides, and which compiler options were used for compiling NEST. If nest-config is not in your path, you need to provided it explicitly like this

```{r, engine='bash'}
cmake -Dwith-nest=${NEST_INSTALL_DIR}/bin/nest-config ../MyModule
```

MyModule will then be installed to `${NEST_INSTALL_DIR}`. This ensures that NEST will be able to find initializing SLI files for the module.

4. Compile.

```{r, engine='bash'}
make
make install
```


### Using MyModule

#### In PyNEST, use

```python
nest.Install("mymodule")
```