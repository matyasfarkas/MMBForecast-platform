##
## This script recursively goes through the matlab directory and its
## subdirectories in the unstable branch, comparing the calling structures of the
## functions to those in another branch/commit, in this case the 4.2.5 branch
## The output is compatible with a Dynare Wiki page
##
## NB: To err on the side of safety, this script won't clean your workspace so a
## git command to checkout a different branch will fail if your workspace is
## dirty.
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
import subprocess
import sys

def runCompare():
    unstable = getFuncDict()
    subprocess.call(['git','checkout','4.2'])

    stable = getFuncDict()
    subprocess.call(['git','checkout','master'])

    us = set(unstable.keys())
    ss = set(stable.keys())

    deletedFunctions = sorted(ss - us)
    makeOldFunctionsWikiTable(deletedFunctions, stable)

    newFunctions = sorted(us - ss)
    makeNewFunctionsWikiTable(newFunctions, unstable)

    outputChanged = set([])
    inputChanged = set([])
    commonFunctions = sorted(ss & us)
    for le in commonFunctions:
        sOut = set(stable[le]['output'])
        uOut = set(unstable[le]['output'])
        sIn = set(stable[le]['input'])
        uIn = set(unstable[le]['input'])

        if len(uOut - sOut):
            outputChanged.add(le)

        if len(uIn - sIn):
            inputChanged.add(le)

    unionChanged = sorted(outputChanged | inputChanged)
    makeUnionChangedWikiTable(unionChanged, stable, unstable)


def makeUnionChangedWikiTable(keys, stable, unstable):
    write = sys.stdout.write
    print "= Functions Whose Arguments Have Changed between Dynare 4.2.5 and " \
        "Dynare 4.3 ="
    print "|| '''Location''' || '''Old Output''' || " \
        "'''New Output''' || '''Old Input''' || '''New Input''' ||"
    for k in keys:
        write('|| {{{' + unstable[k]['filename'][10:]  + '}}} || ')
        write(str(stable[k]['output']) + ' || ')
        write(str(unstable[k]['output']) + ' || ')
        write(str(stable[k]['input']) + ' || ')
        write(str(unstable[k]['input']) + ' ||')
        print


def makeNewFunctionsWikiTable(keys, dictionary):
    print '= New Functions in Dynare 4.3 ='
    makeWikiTable(keys, dictionary)


def makeOldFunctionsWikiTable(keys, dictionary):
    print '= Functions Removed in Dynare 4.3 ='
    makeWikiTable(keys, dictionary)


def makeWikiTable(keys, dictionary):
    write = sys.stdout.write
    print "|| '''Location''' || '''Output''' || '''Input''' ||"
    for k in keys:
        write('|| {{{' + dictionary[k]['filename'][10:]  + '}}} || ')
        write(str(dictionary[k]['output']) + ' || ')
        write(str(dictionary[k]['input']) + ' ||')
        print


def printDict(dictionary, title):
    print
    print '***********************************'
    print '** ' + title
    print '***********************************'
    for key in sorted(dictionary.iterkeys()):
        b = dictionary[key]
        print key + ':' + b['filename']
    print '***********************************'
    print '***********************************'
    print


def getFuncDict():
    functions = {}

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
                if inComment:
                    if funcDef[-1:] == '\n':
                        inComment = False
                else:
                    if funcDef[-1:] == '%':
                        inComment = True
                    elif funcDef[-1:] == ')':
                        break
            funcDef = funcDef.strip()
            f.close()

            funcDict = {}
            funcDict['filename'] = fullfilename

            # OUTPUT
            spliteq = string.rsplit(funcDef, '=')
            if len(spliteq) == 1:
                # no output
                outputs = ['']
                funcDict['output'] = outputs
                rhs = string.split(funcDef, ' ', 1).pop(1).strip()
            else:
                outputs = string.split(spliteq.pop(0), ' ', 1).pop(1)
                outputs = [string.strip(outputs, '[] ')]
                outputs = getputs(outputs)
                funcDict['output'] = outputs
                rhs = string.strip(spliteq.pop(0), ' .\n')

            # FUNCTION NAME
            splitfn = string.split(rhs, '(')
            fn = splitfn.pop(0)

            # INPUT
            inputs = [string.strip(le,')') for le in splitfn]
            inputs = getputs(inputs)
            funcDict['input'] = inputs

            functions[fn] = funcDict

    return functions


def getputs(puts):
    ''' Can be of the forms
    a
    a
    a, b
    a, b   % aoeu
    a, ...
       b
    a, ... % aoeu
       b

    Return a list of (in/out puts)
    '''
    splitout = string.split(string.join(puts), '\n')

    for i,le in enumerate(splitout[:]):
        indx = string.find(le, '%')
        if indx != -1:
            splitout[i] = string.strip(le[:indx], '\t .')
        else:
            splitout[i] = string.strip(le, '\t .')

    splitcommas = string.split(string.join(splitout), ',')
    return [string.strip(le) for le in splitcommas]


if __name__ == "__main__":
    runCompare()
