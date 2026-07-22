# **************************************************************************
# *
# * Authors:   Blanca Pueche (blanca.pueche@cnb.csis.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import json
import shutil

import os, json
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pwem.protocols import EMProtocol
from pyworkflow.object import String, Float

from openstructure.__init__ import Plugin
from pwchem.__init__ import Plugin as pwchemPlugin
from pwem.objects.data import AtomStruct, SetOfAtomStructs
from pwchem.objects.base import SmallMolecule, SetOfSmallMolecules
from pwem.convert import cifToPdb

from openstructure import OPENSTRUCT_DIC
from pwchem.constants import OPENBABEL_DIC


class LigandSelect(Select):
    def accept_residue(self, residue):
        # Keep only hetero residues (ligands), skip waters
        return residue.id[0].startswith("H_") and residue.get_resname() != "HOH"


class ProtCompareStructures(EMProtocol):
    """
    AI Generated:

    This protocol compares a predicted macromolecular structure against a
    reference structure using OpenStructure (OST). It computes a comprehensive
    set of global and interface-level structural similarity metrics, generates
    superposed structures, and annotates the predicted model with per-residue
    Local LDDT scores for visualization and downstream analysis.

    The protocol supports both monomeric and multimeric structures and can
    optionally compute backbone-only LDDT and quaternary structure similarity
    (QS-score).

    Core Concepts
    -------------
    LDDT:
        Local Distance Difference Test, a superposition-independent metric that
        measures local structural agreement between the predicted and reference
        structures.

    Local LDDT:
        Per-residue LDDT score describing the local accuracy of each residue.
        The protocol stores these values in the occupancy column of the output
        PDB file, allowing direct visualization in molecular graphics software.

    TM-score:
        Length-independent measure of global structural similarity between two
        structures.

    DockQ:
        Interface quality metric for protein complexes, combining interface
        RMSD, ligand RMSD, and native contact recovery into a single score.

    QS-score:
        Quaternary Structure Score measuring the similarity between the
        interfaces and overall organization of two macromolecular complexes.

    GDT Scores:
        Global Distance Test metrics (GDT-TS and GDT-HA) evaluating structural
        similarity under different distance thresholds.

    Workflow
    --------
    1. Convert input structures to PDB format if necessary.
    2. Perform structural comparison using OpenStructure.
    3. Compute global structural similarity metrics.
    4. Optionally compute backbone-only LDDT.
    5. Optionally compute QS-score for multimeric assemblies.
    6. Compute per-residue Local LDDT values.
    7. Write Local LDDT scores into the occupancy field of the predicted model.
    8. Export the superposed model and reference structures.

    Input
    -----
    - inputModel:
        Predicted structure to evaluate.

    - inputReference:
        Reference (native or experimental) structure.

    Parameters
    ----------
    - Chain mapping:
        Optional correspondence between chains of the predicted and reference
        structures. Useful when chain identifiers differ or when only selected
        interfaces should be evaluated.

    - Generate backbone LDDT:
        Computes LDDT using only backbone atoms (CA for proteins and C3' for
        nucleic acids).

    - Compute QS-score:
        Calculates the Quaternary Structure Score to evaluate similarity
        between multimeric assemblies.

    Computed Metrics
    ----------------
    Depending on the selected options, the protocol reports:

    - Global LDDT
    - Backbone LDDT
    - TM-score
    - DockQ (average and weighted)
    - QS-score
    - GDT-TS
    - GDT-HA
    - Global RMSD after rigid-body superposition

    For multimeric complexes, the protocol also reports per-interface:

    - QS-score
    - DockQ
    - Interface RMSD (iRMSD)

    Output
    ------
    - outputAtomStruct:
        Predicted structure after superposition onto the reference.

    The model contains:
        - Per-residue Local LDDT stored in the occupancy column.
        - Global comparison metrics (LDDT, TM-score, RMSD, DockQ, QS-score, etc.)
            stored as object attributes for downstream protocols.

    The predicted model contains the Local LDDT score of every residue stored
    in the occupancy field, enabling residue-level coloring in molecular
    visualization software such as ChimeraX or PyMOL.

    Additional Files
    ----------------
    - compare_structures.json:
        JSON report containing all global, interface-level, and per-residue
        comparison metrics produced by OpenStructure.

    Use Cases
    ---------
    - Evaluating predicted structures against experimental references
    - Benchmarking protein structure prediction methods
    - Assessing protein complex and multimer accuracy
    - Comparing refined models
    - Visualizing residue-level prediction accuracy using Local LDDT
    - Measuring interface similarity in macromolecular assemblies
    """
    _label = 'compare structures'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputModel', params.PointerParam, allowsNull=True,
                      pointerClass='AtomStruct',
                      label="Input predicted model: ",
                      help='Select the predicted model.')

        form.addParam('inputReference', params.PointerParam, allowsNull=False,
                      pointerClass='AtomStruct',
                      label="Input reference structure: ",
                      help='Select the reference structure.')


        group = form.addGroup('Parameters')
        group.addParam('mapping', params.StringParam, default='', expertLevel = LEVEL_ADVANCED,
                       label='Chain mapping',
                       help='Specify a chain mapping between model and native structure. If the native contains two chains "H" and "L" while the model contains two chains "A" and "B", and chain A is a model of native chain H and chain B is a model\n '
                            'of native chain L, the flag can be set as: "--mapping AB:HL". This can also help limit the search to specific native interfaces. For example, if the native is a tetramer (ABCD) but the user is only interested in the \n'
                            'interface between chains B and C, the flag can be set as: "--mapping :BC" or the equivalent "--mapping *:BC".'
                       )
        group.addParam('backboneLDDT', params.BooleanParam, default=True,
                       label='Generate backbone lddt: ',
                       help="LDDT in this case is only computed on backbone atoms: CA for peptides and C3' for nucleotides."
                       )
        group.addParam('qsScore', params.BooleanParam, default=True,
                       label='Compute QS score: ',
                       help="The QS-score (Quaternary Structure Score) is a metric designed to compare the quaternary structure of two macromolecular complexes. \n"
                            "Evaluates how similar the interfaces between chains are."
                       )

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertFilesStep)
        self._insertFunctionStep(self.runOSTProtStep)
        self._insertFunctionStep(self.writeLocalLDDTToModelStep)

        self._insertFunctionStep(self.createOutputStep)

    def convertFilesStep(self):
        inModel = self.inputModel.get().getFileName()
        inpPDBModel = self._getExtraPath("model.pdb")
        self.convertOrCopy(inModel, inpPDBModel)
        inRef = self.inputReference.get().getFileName()
        inpPDBRef = self._getExtraPath("reference.pdb")
        self.convertOrCopy(inRef, inpPDBRef)

    def runOSTProtStep(self):
        args = []
        args.extend([
            "-m", str(os.path.abspath(self._getExtraPath("model.pdb"))),
            "-r", str(os.path.abspath(self._getExtraPath("reference.pdb"))),
            "-o", str(os.path.abspath(self._getPath('compare_structures.json'))),
            '-d',
            '--lddt',
            '--dockq',
            '--tm-score',
            '--rigid-scores', #GDT scores
            '--local-lddt'
        ])

        if self.mapping.get() != '':
            args += ['-c', self.mapping.get()]
        if self.backboneLDDT.get():
            args.append('--bb-lddt')
        if self.qsScore.get():
            args.append('--qs-score')

        Plugin.runCondaCommand(
            self,
            args=" ".join(args),
            condaDic=OPENSTRUCT_DIC,
            program="ost compare-structures",
            cwd=os.path.abspath(Plugin.getVar(OPENSTRUCT_DIC['home']))
        )

    def writeLocalLDDTToModelStep(self):
        jsonFile = self._getPath("compare_structures.json")
        modelFile = self._getExtraPath("model_compare_structures.pdb")

        if not os.path.exists(jsonFile):
            raise FileNotFoundError(jsonFile)

        if not os.path.exists(modelFile):
            raise FileNotFoundError(modelFile)

        with open(jsonFile) as f:
            data = json.load(f)

        local_lddt = data.get("local_lddt", {})

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("model", modelFile)

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] != " ":
                        continue
                    chainId = chain.id
                    resNum = residue.id[1]
                    insCode = residue.id[2].strip()

                    key = f"{chainId}.{resNum}.{insCode}"
                    score = local_lddt.get(key)
                    if score is None:
                        score = 0.0
                    for atom in residue:
                        atom.set_occupancy(float(score))
        io = PDBIO()
        io.set_structure(structure)
        io.save(modelFile)

    def createOutputStep(self):
        jsonFile = self._getPath("compare_structures.json")
        stats = {}
        if os.path.exists(jsonFile):
            with open(jsonFile) as f:
                stats = json.load(f)

        modelFile = self._getExtraPath("model_compare_structures.pdb")
        if not os.path.exists(modelFile):
            raise FileNotFoundError(modelFile)

        model = self.inputModel.get().clone()
        model.setFileName(modelFile)

        for jsonKey, attrName in [
            ("lddt", "_lddt"),
            ("bb_lddt", "_bbLddt"),
            ("tm_score", "_tmScore"),
            ("qs_global", "_qsScore"),
            ("dockq_ave", "_dockQAve"),
            ("dockq_wave", "_dockQWeighted"),
            ("oligo_gdtts", "_gdtTs"),
            ("oligo_gdtha", "_gdtHa"),
            ("rmsd", "_rmsd"),
        ]:
            value = stats.get(jsonKey)
            if value is not None:
                setattr(model, attrName, Float(value))

        self._defineOutputs(outputAtomStruct=model)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []

        json_path = os.path.abspath(self._getPath('compare_structures.json'))

        if not os.path.exists(json_path):
            return ["Comparison JSON file not found"]

        with open(json_path, "r") as f:
            results = json.load(f)

        def format_value(key):
            value = results.get(key)
            if isinstance(value, (int, float)):
                return f"{value:.3f}"
            return str(value) if value is not None else "N/A"

        summary.append(f"LDDT: {format_value('lddt')}")
        summary.append(f"Backbone LDDT: {format_value('bb_lddt')}")
        summary.append(f"TM-score: {format_value('tm_score')}")
        summary.append(f"QS-score (global): {format_value('qs_global')}")
        summary.append(f"DockQ (average): {format_value('dockq_ave')}")
        summary.append(f"DockQ (weighted): {format_value('dockq_wave')}")
        summary.append(f"GDT-TS: {format_value('oligo_gdtts')}")
        summary.append(f"GDT-HA: {format_value('oligo_gdtha')}")
        summary.append(f"Global RMSD (after rigid superposition): {format_value('rmsd')} Å")

        chains = results.get("model_chains", [])
        summary.append(f"Number of chains: {len(chains)}")

        qs_interfaces = results.get("qs_interfaces", [])
        qs_scores = results.get("per_interface_qs_global", [])
        dockq_scores = results.get("dockq", [])
        irmsd_scores = results.get("irmsd", [])

        if qs_interfaces:
            summary.append("\nPer-interface QS scores:")
            for interface, qs in zip(qs_interfaces, qs_scores):
                summary.append(
                    f"{interface[0]}-{interface[1]}: QS={qs:.3f}"
                )

        dockq_interfaces = results.get("dockq_interfaces", [])
        dockq_scores = results.get("dockq", [])
        irmsd_scores = results.get("irmsd", [])

        if dockq_interfaces:
            summary.append("\nPer-interface DockQ scores:")
            for interface, dockq, irmsd in zip(
                    dockq_interfaces, dockq_scores, irmsd_scores):
                summary.append(
                    f"{interface[0]}-{interface[1]}: "
                    f"DockQ={dockq:.3f},    iRMSD={irmsd:.3f} Å"
                )

        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        #dockedProtein = os.path.splitext(
        #    os.path.basename(os.path.abspath(self.inputLigands.get().getProteinFile()))
        #)[0]
        #inputModel = os.path.splitext(
        #    os.path.basename(os.path.abspath(self.inputModel.get().getFileName()))
        #)[0]
        validations = []
        #if dockedProtein != inputModel:
        #    validations.append('Ligand selected is not docked in model structure.')
        return validations

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def convertOrCopy(self, inModel, inpPDBModel):
        if inModel.endswith('.cif'):
            cifToPdb(inModel, inpPDBModel)
        else:
            shutil.copy(inModel, inpPDBModel)

    def getDockedLigand(self):
        myMol = None
        for mol in self.inputLigands.get():
            if mol.__str__() == self.inputMolecule.get():
                myMol = mol.clone()
                break

        if myMol is None:
            print('The input ligand is not found')
        return myMol

    def extractLigand(self, refFile, outFile):
        if refFile.endswith(".cif"):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        structure = parser.get_structure("ref", refFile)

        io = PDBIO()
        io.set_structure(structure)
        io.save(outFile, LigandSelect())
