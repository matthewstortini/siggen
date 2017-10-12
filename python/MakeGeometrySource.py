#!/usr/local/bin/python
from __future__ import print_function
from shutil import copyfile
import fileinput

#Its impossible to wrap a c++ template class in a cython object... so this manually copy-pastes source code for each geometry type

def make_source(geometry_list):

    source_file = "siggen/PySiggen.pyi.src"

    for geom in geometry_list:
        new_file = "siggen/PySiggen_%s.pxi" % geom
        copyfile(source_file, new_file)

        f = fileinput.FileInput(new_file, inplace=True)
        for line in f:
            print(line.replace("{DetectorType}", geom), end='')
        f.close()
