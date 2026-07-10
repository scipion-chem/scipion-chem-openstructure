import os

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb, ProtImportVolumes

from .. import Plugin
from ..protocols import ProtEMProt
from ..utils import assertHandle

class TestEMProt(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls._runImportPDB()
        cls._runImportVolume()

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0,
            pdbId='9r1h',
        )
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runImportVolume(cls):
        protImpVol = cls.newProtocol(
            ProtImportVolumes,
            importFrom=1,
            emdbId='53509'
        )
        cls.launchProtocol(protImpVol)
        cls.protImpVol = protImpVol

    def _runEMProt(self):
        protEMProt = self.newProtocol(ProtEMProt)
        protEMProt.useGpu.set('False')

        protEMProt.inputVolume.set(self.protImpVol.outputVolume)
        protEMProt.inputStructure.set(self.protImportPDB.outputPdb)

        self.proj.launchProtocol(protEMProt, wait=True)
        return protEMProt

    def test(self):
        protEMProt = self._runEMProt()
        self._waitOutput(protEMProt, 'outputAtomStruct', sleepTime=5)

        assertHandle(self.assertIsNotNone, getattr(protEMProt, 'outputAtomStruct', None), cwd=protEMProt.getWorkingDir())
