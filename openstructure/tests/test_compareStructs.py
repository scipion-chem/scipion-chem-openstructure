import os

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportSetOfAtomStructs

from .. import Plugin
from ..protocols import ProtCompareStructures
from ..utils import assertHandle

class TestCompareStructures(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls._runImportPDB()
        cls._waitOutput(cls.protImportPDB, 'outputAtomStructs', sleepTime=5)

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportSetOfAtomStructs,
            inputPdbData=0,
            pdbIds='9j42, 9j4a',
        )
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB


    def _runOST(self):
        protOST = self.newProtocol(ProtCompareStructures)

        protOST.inputModel.set(self.protImportPDB.outputAtomStructs[0])
        protOST.inputReference.set(self.protImportPDB.outputAtomStructs[1])

        self.proj.launchProtocol(protOST, wait=True)
        return protOST

    def test(self):
        protOST = self._runOST()
        self._waitOutput(protOST, 'outputAtomStructs', sleepTime=5)

        assertHandle(self.assertIsNotNone, getattr(protOST, 'outputAtomStructs', None), cwd=protOST.getWorkingDir())
