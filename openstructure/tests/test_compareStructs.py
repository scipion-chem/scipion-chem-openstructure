import os

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb

from .. import Plugin
from ..protocols import ProtCompareStructures
from ..utils import assertHandle

class TestCompareStructures(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls._runImportPDBModel()
        cls._waitOutput(cls.protImportPDB1, 'outputPdb', sleepTime=5)
        cls._runImportPDBRef()
        cls._waitOutput(cls.protImportPDB2, 'outputPdb', sleepTime=5)

    @classmethod
    def _runImportPDBModel(cls):
        protImportPDB1 = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0,
            pdbId='9j42',
        )
        cls.launchProtocol(protImportPDB1)
        cls.protImportPDB1 = protImportPDB1

    @classmethod
    def _runImportPDBRef(cls):
        protImportPDB2 = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0,
            pdbId='9j4a',
        )
        cls.launchProtocol(protImportPDB2)
        cls.protImportPDB2 = protImportPDB2


    def _runOST(self):
        protOST = self.newProtocol(ProtCompareStructures)

        protOST.inputModel.set(self.protImportPDB1.outputPdb)
        protOST.inputReference.set(self.protImportPDB2.outputPdb)

        self.proj.launchProtocol(protOST, wait=True)
        return protOST

    def test(self):
        protOST = self._runOST()
        self._waitOutput(protOST, 'outputAtomStruct', sleepTime=5)

        assertHandle(self.assertIsNotNone, getattr(protOST, 'outputAtomStruct', None), cwd=protOST.getWorkingDir())
