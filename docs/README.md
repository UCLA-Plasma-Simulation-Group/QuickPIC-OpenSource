# Generate Documentation

QuickPIC uses [FORD](https://github.com/cmacmackin/ford) to generate documentation. Before installing FORD, make sure Python and relevant libraries are correctly installed. The simplest way to install FORD is using [pip](https://pip.pypa.io/en/latest/) with this command

```
sudo pip install ford
```

which will install FORD and all its [dependencies](https://github.com/cmacmackin/ford/wiki/Dependencies) automatically. If you wish to install the most recent version (in case that FORD installed through pip runs abnormally in some systems), you can also clone the github repository:

```
git clone https://github.com/cmacmackin/ford.git
```

or download the [install package](https://github.com/cmacmackin/ford/archive/master.zip). Go into the FORD folder and run the command below to install manually:

```
python setup.py install
```

If you wish to generate relation graphs in the documentation, the Graphviz binary is needed. On Debian based systems, this can be done with

```
sudo apt-get install graphviz
```

The project file (quickpic_ford.md) includes various [options](https://github.com/cmacmackin/ford/wiki/Project-File-Options) and information of the QuickPIC documentation in the meta-data. Be aware that the option __include__, which defines the directories for file searching of Fortran's intrinsic `include` statement, should be properly configured before generating the documentation. In current version of QuickPIC, only the path of `mpif.h` needs to be added to __include__ option in FORD_DOC.md, just as below

```
include: [path of "mpif.h"]
```

Before using FORD, keep in mind that you must make sure the QuickPIC source code compiles correctly because FORD is not made to check any syntax error. Run the command below to generate the documentation:

```
ford quickpic_ford.md
```

# Online documentation

Visit [here](https://ucla-plasma-simulation-group.github.io/QuickPIC-OpenSource/) to see the documentation.