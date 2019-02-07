##
## This script recursively goes through the matlab directory and its
## subdirectories, comparing all .m file names to the function declarations at
## the top of those files. If the two do not correspond, it prints out the name
## of the file and the name of the function declared in the file.
##

##
## Copyright (C) 2012 Dynare Team
##
## This file is part of Dynare.
##
## Dynare is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Dynare is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
##

import os
import string

for dirname, dirnames, filenames in os.walk('../matlab'):
    for filename in filenames:
        filename = string.strip(filename)

        if filename[-2:] != '.m' or filename == 'msstart2.m' or filename == 'msstart_setup.m' or filename == 'qmc_sequence.m':
            continue

        fullfilename = os.path.join(dirname, filename)
        f = open(fullfilename, 'r')
        funcDef = ''
        inComment = False
        while True:
            funcDef += f.read(1)
            if funcDef[-1:] == '%':
                inComment = True

            if inComment:
                if funcDef[-1:] == '\n':
                    inComment = False
            else:
                if funcDef[-1:] == '(':
                    break
        f.close()

        spliteq = string.rsplit(funcDef, '=')
        if len(spliteq) == 1:
            spliteq = string.rsplit(funcDef, 'function ')

        spliteq = spliteq.pop()
        spliteq = string.strip(spliteq, '. ')
        spliteq = string.strip(spliteq, '\n ')
        spliteq = string.strip(spliteq, '( ')

        if filename[:-2] != spliteq:
            print fullfilename + ': ' + spliteq

