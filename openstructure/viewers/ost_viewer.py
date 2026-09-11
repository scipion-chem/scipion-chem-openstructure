from tkinter import messagebox
import os

from pwem.viewers import Chimera
from pyworkflow.viewer import Viewer
import pyworkflow.viewer as pwviewer
from pwchem.viewers import PyMolViewer

from openstructure.protocols import ProtCompareStructures
from Bio.PDB import PDBParser

from pwchem.__init__ import Plugin as pwchemPlugin

from pyworkflow.protocol.params import LabelParam, BooleanParam, EnumParam
from pyworkflow.viewer import Viewer
from pwem.viewers import Chimera
import json
import tempfile
import matplotlib.pyplot as plt


class ProtCompareStructuresViewer(pwviewer.ProtocolViewer):

    _label = "viewer compare structures"
    _targets = [ProtCompareStructures]


    def __init__(self, **args):
        super().__init__(**args)
        self.CHOICES = self._getInterfaceChoices()


    def _defineParams(self, form):

        form.addSection(label='Visualization')
        form.addParam('showAlignment',
                      LabelParam,
                      label='Open aligned structures ')
        form.addParam('colorLDDT',
                      BooleanParam, default=True,
                      label='Color by local lDDT')

        group = form.addGroup('Graphs')
        group.addParam('showLocalLDDT',
                      LabelParam,
                      label='Plot local lDDT')
        group.addParam('showInterfaces',
                      LabelParam,
                      label='Interface scores')
        group.addParam('showSummary',
                      LabelParam,
                      label='Summary table')

        group = form.addGroup('Visualize with PLIP')
        group.addParam('interfaceChoice', EnumParam, default=0, choices=self.CHOICES,
                       label='Display chain interactions: ',
                       help='Display this interface interactions.')
        form.addParam('showPLIPinterface',
                      LabelParam,
                      label='Open PLIP ')

    def _getVisualizeDict(self):
        return {
            'showAlignment': self._viewAlignment,
            'showLocalLDDT': self._viewLocalLDDT,
            'showInterfaces': self._viewInterfaces,
            'showSummary': self._viewSummary,
            'showPLIPinterface': self._viewPLIPInterface,
        }

    def _getInterfaceChoices(self):
        jsonFile = self.protocol._getPath("compare_structures.json")
        if not os.path.exists(jsonFile):
            return ["No interfaces"]
        with open(jsonFile) as f:
            data = json.load(f)
        chains = data.get("model_chains", [])
        return chains if chains else ["No interfaces"]

    def _viewAlignment(self, paramName=None):
        if self.colorLDDT.get():
            return self._viewAlignmentLDDT()
        else:
            return self._viewAlignmentPlain()

    def _viewAlignmentPlain(self, paramName=None):
        fnCmd = self.protocol._getExtraPath("chimera_compare.cxc")

        modelFile = os.path.abspath(
            self.protocol.outputAtomStruct.getFileName()
        )
        refFile = os.path.abspath(
            self.protocol._getExtraPath("reference_compare_structures.pdb")
        )

        with open(fnCmd, "w") as f:
            f.write(f"open {refFile}\n")
            f.write(f"open {modelFile}\n")

            f.write("hide atoms\n")
            f.write("cartoon\n")
            f.write("transparency #1 60\n")
            f.write("match #2 to #1\n")
            f.write("view orient\n")

        return [Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")]

    def _viewAlignmentLDDT(self, paramName=None):

        fnCmd = self.protocol._getExtraPath("chimera_lddt.cxc")

        modelFile = os.path.abspath(
            self.protocol.outputAtomStruct.getFileName()
        )
        refFile = os.path.abspath(
            self.protocol._getExtraPath("reference_compare_structures.pdb")
        )

        with open(fnCmd, "w") as f:
            f.write(f"open {refFile}\n")
            f.write(f"open {modelFile}\n")

            f.write("hide atoms\n")
            f.write("cartoon\n")

            f.write("color #1 gray\n")
            f.write("transparency #1 60\n")

            f.write("match #2 to #1\n")
            f.write("color byattribute occupancy palette white:mintcream:lightgreen:forestgreen:darkgreen\n")
            f.write("key white:0 lightgreen:0.5 darkgreen:1\n")
            f.write("color #1 cornflowerblue\n")
            f.write("view orient\n")

        return [Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")]

    def _viewLocalLDDT(self, param=None):

        jsonFile = self.protocol._getPath("compare_structures.json")

        with open(jsonFile) as f:
            data = json.load(f)

        local = data["local_lddt"]

        chains = {}

        for key, value in local.items():
            if value is None:
                continue

            chain, res, _ = key.split('.')

            chains.setdefault(chain, {"x": [], "y": []})

            chains[chain]["x"].append(int(res))
            chains[chain]["y"].append(value)

        plt.figure(figsize=(10, 4))

        for chain in sorted(chains):
            plt.plot(chains[chain]["x"],
                     chains[chain]["y"],
                     label=f"Chain {chain}")

        plt.xlabel("Residue")
        plt.ylabel("Local lDDT")
        plt.ylim(0, 1)
        plt.grid(True)
        plt.legend()

        plt.show()

    def _viewInterfaces(self, param=None):
        jsonFile = self.protocol._getPath("compare_structures.json")
        with open(jsonFile) as f:
            data = json.load(f)

        from tkinter import messagebox
        import numpy as np

        qs_interfaces = data.get("qs_interfaces", [])
        dockq_interfaces = data.get("dockq_interfaces", [])

        if not qs_interfaces and not dockq_interfaces:
            messagebox.showwarning(
                "No interface scores",
                "No interfaces were detected.\n\n"
                "QS-score, DockQ and iRMSD are only available for multimeric complexes."
            )
            return []

        qs = data.get("per_interface_qs_global", [])
        dockq = data.get("dockq", [])
        irmsd = data.get("irmsd", [])

        fig, (ax1, ax2) = plt.subplots(
            2, 1,
            figsize=(max(10, len(qs_interfaces) * 0.4), 8)
        )

        # ------------------------
        # QS-score
        # ------------------------
        if qs_interfaces:
            labels = [f"{i[0]}-{i[1]}" for i in qs_interfaces]
            x = np.arange(len(labels))

            ax1.bar(x, qs)
            ax1.set_title("Per-interface QS-score")
            ax1.set_ylabel("QS")
            ax1.set_xticks(x)
            ax1.set_xticklabels(labels, rotation=60, ha="right")
            ax1.set_ylim(0, 1.05)

        # ------------------------
        # DockQ + iRMSD
        # ------------------------
        if dockq_interfaces:
            labels = [f"{i[0]}-{i[1]}" for i in dockq_interfaces]
            x = np.arange(len(labels))
            w = 0.35

            ax2.bar(x - w / 2, dockq, width=w, label="DockQ")
            ax2.bar(x + w / 2, irmsd, width=w, label="iRMSD (Å)")

            ax2.set_title("Per-interface DockQ / iRMSD")
            ax2.set_ylabel("Score")
            ax2.set_xticks(x)
            ax2.set_xticklabels(labels, rotation=60, ha="right")
            ax2.legend()

        plt.tight_layout()
        plt.show()

    def _viewSummary(self, param=None):

        jsonFile = self.protocol._getPath("compare_structures.json")

        with open(jsonFile) as f:
            data = json.load(f)

        metrics = [
            ("lDDT", data.get("lddt")),
            ("BB-lDDT", data.get("bb_lddt")),
            ("TM-score", data.get("tm_score")),
            ("QS-score", data.get("qs_global")),
            ("DockQ", data.get("dockq_ave")),
            ("GDT-TS", data.get("oligo_gdtts")),
            ("GDT-HA", data.get("oligo_gdtha"))
        ]

        metrics = [(n, v) for n, v in metrics if v is not None]

        names = [m[0] for m in metrics]
        values = [m[1] for m in metrics]

        plt.figure(figsize=(8, 4))
        plt.barh(names, values)

        plt.xlim(0, 1)
        plt.xlabel("Score")
        plt.title("Global structure comparison")

        for i, v in enumerate(values):
            plt.text(v + 0.02, i, f"{v:.3f}")

        plt.tight_layout()
        plt.show()

    def _viewPLIPInterface(self, param=None):
        chain = self._getInterfaceChoices()[self.interfaceChoice.get()]

        stype = self._getStructureType()
        pmlDir = os.path.abspath(self.protocol._getExtraPath("plip_interface"))
        os.makedirs(pmlDir, exist_ok=True)

        modelFile = os.path.abspath(
            self.protocol.outputAtomStruct.getFileName()
        )

        if stype == "dna" or stype == "rna":
            args = f"-f {modelFile} --dnareceptor --inter {chain} -y -o {pmlDir}"
        else:
            args = f"-f {modelFile} --inter {chain} -y -o {pmlDir}"

        pwchemPlugin.runPLIP(args, cwd=pmlDir)

        pseFile = None
        for root, _, files in os.walk(pmlDir):
            for f in files:
                if f.endswith(".pse"):
                    pseFile = os.path.join(root, f)
                    break
            if pseFile:
                break


        if pseFile is None:
            messagebox.showwarning(
                "PLIP",
                f"No interactions were found for {chain}."
            )
            return []

        pymolV = PyMolViewer(project=self.getProject())
        return pymolV._visualize(
            os.path.abspath(pseFile),
            cwd=os.path.dirname(pseFile)
        )

    def _getStructureType(self):
        modelFile = os.path.abspath(
            self.protocol.outputAtomStruct.getFileName()
        )

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("model", modelFile)

        hasDNA = False
        hasRNA = False

        DNA_RES = {"DA", "DT", "DG", "DC", "DI"}
        RNA_RES = {"A", "U", "G", "C", "I"}

        for residue in structure.get_residues():
            if residue.id[0] != " ":
                continue

            name = residue.resname.strip()

            if name in DNA_RES:
                hasDNA = True
            elif name in RNA_RES:
                hasRNA = True

        if hasDNA:
            return "dna"
        if hasRNA:
            return "rna"
        return "protein"