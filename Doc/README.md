These documentions may be build with the [`arara`](https://github.com/cereda/arara) automated build system. `arara` is part of [`CTAN`](https://ctan.org/pkg/arara). To compile the documents use

    arara [PROJECTNAME].tex

e.g.

    arara Peridigm_Users_Guide.tex

in the respective document main path.

This type of compilation is currently only supported under Windows due to the required batch files. In case of problems, please also note the information in the respective headers of the project files.

Each of the individual documentations may be build seperately. However, all depend on definitions from the `PeriDoX_Common` folder. Thus, you need this folder in the exact relative position to where it can be found in this repository.
