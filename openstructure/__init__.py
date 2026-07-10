# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os
from os.path import join, exists

import pwem
from scipion.install.funcs import InstallHelper

from pyworkflow import SPA, TOMO, MODELLING
from .constants import OPENSTRUCT_DIC

_version_ = '2.11'
_logo = "logo.png"
_references = ['']




class Plugin(pwem.Plugin):
    _homeVar = OPENSTRUCT_DIC['home']
    _pathVars = [OPENSTRUCT_DIC['home']]
    _supportedVersions = [OPENSTRUCT_DIC['version']]

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(OPENSTRUCT_DIC['home'], OPENSTRUCT_DIC['name'] + '-' + OPENSTRUCT_DIC['version'])

    @classmethod
    def defineBinaries(cls, env, default=True):

        installer = InstallHelper(
            OPENSTRUCT_DIC['name'],
            packageHome=cls.getVar(OPENSTRUCT_DIC['home']),
            packageVersion=OPENSTRUCT_DIC['version']
        )

        installer.getCondaEnvCommand(
            OPENSTRUCT_DIC['name'],
            binaryVersion=OPENSTRUCT_DIC['version'],
            pythonVersion='3.11'
        ).addCommand(
            f"{cls.getEnvActivationCommand(OPENSTRUCT_DIC)} && "
            "conda install -y -c conda-forge -c bioconda openstructure",
            "OPENSTRUCTURE_INSTALLED"
        ).addPackage(
            env,
            dependencies=['git', 'conda'],
            default=default
        )

    @classmethod
    def getEnvName(cls, packageDictionary):
        """ This function returns the name of the conda enviroment for a given package. """
        return '{}-{}'.format(packageDictionary['name'], packageDictionary['version'])

    @classmethod
    def getEnvActivationCommand(cls, packageDictionary, condaHook=True):
        """ This function returns the conda enviroment activation command for a given package. """
        return '{}conda activate {}'.format(cls.getCondaActivationCmd() if condaHook else '',
                                            cls.getEnvName(packageDictionary))

    @classmethod
    def runCondaCommand(cls, protocol, args, condaDic, program, cwd=None, popen=False, silent=True, retOut=False):
        """ General function to run conda commands """
        result = None
        fullProgram = f'{cls.getEnvActivationCommand(condaDic)} && {program} '
        if not popen and not retOut:
            protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd, numberOfThreads=1)
        else:
            if not retOut:
                kwargs = {}
                if silent:
                    kwargs = {"stdout": subprocess.DEVNULL, "stderr": subprocess.DEVNULL}
                run(fullProgram + args, env=cls.getEnviron(), cwd=cwd, shell=True, **kwargs)
            else:
                result = subprocess.check_output(fullProgram + args, cwd=cwd, shell=True, text=True)
        return result


