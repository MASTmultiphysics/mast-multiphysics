MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
Copyright (C) 2013-2020  Manav Bhatia and MAST authors

[![Build Status](https://travis-ci.com/MASTmultiphysics/mast-multiphysics.svg?branch=master)](https://travis-ci.com/MASTmultiphysics/mast-multiphysics)

This code was developed under funding from the Air Force Research Laboratory. 
MAST was cleared for public release on 08 Nov 2016 with case number 88ABW-2016-5689. 

Documentation for the code is available at [https://mastmultiphysics.github.io](https://mastmultiphysics.github.io).

## Submodules
To keep the size of this main MAST repository smaller, a git submodule is used 
to store larger media/assets such as images and animations used for the documentation 
in a separate repo (doc/assets). To build the documentation locally, you must update
the submodule. To do this, simply run the following commands from inside the root
level of this main repository:
```
git submodule init
git submodule update
```
