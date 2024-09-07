import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget
from PyQt5.QtCore import Qt, QRect
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


class Molecule3DViewer(QMainWindow):
    def __init__(self, smiles_code):
        super().__init__()

        # Set up the main window
        self.setWindowTitle("3D Molecule Viewer")
        self.setGeometry(100, 100, 800, 600)
        self.snap_to_right()

        # Generate the molecule from SMILES
        self.molecule = Chem.MolFromSmiles(smiles_code)
        if self.molecule is None:
            raise ValueError("Invalid SMILES code")

        AllChem.EmbedMolecule(self.molecule)
        AllChem.UFFOptimizeMolecule(self.molecule)

        # Set up the VTK renderer and interactor
        self.vtk_widget = QVTKRenderWindowInteractor(self)
        self.layout = QVBoxLayout()
        self.layout.addWidget(self.vtk_widget)

        central_widget = QWidget(self)
        central_widget.setLayout(self.layout)
        self.setCentralWidget(central_widget)

        self.renderer = vtk.vtkRenderer()
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)
        self.iren = self.vtk_widget.GetRenderWindow().GetInteractor()

        self.setup_molecule()

        self.show()
        self.iren.Initialize()

    def setup_molecule(self):
        # Convert RDKit molecule to VTK objects
        conformer = self.molecule.GetConformer()
        atom_positions = []

        for i in range(self.molecule.GetNumAtoms()):
            pos = conformer.GetAtomPosition(i)
            atom_positions.append([pos.x, pos.y, pos.z])

        # Render atoms
        for pos in atom_positions:
            sphere = vtk.vtkSphereSource()
            sphere.SetRadius(0.4)
            sphere.SetCenter(pos)
            sphere_mapper = vtk.vtkPolyDataMapper()
            sphere_mapper.SetInputConnection(sphere.GetOutputPort())
            sphere_actor = vtk.vtkActor()
            sphere_actor.SetMapper(sphere_mapper)
            self.renderer.AddActor(sphere_actor)

        # Render bonds
        for bond in self.molecule.GetBonds():
            begin_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            pos1 = atom_positions[begin_atom_idx]
            pos2 = atom_positions[end_atom_idx]

            line = vtk.vtkLineSource()
            line.SetPoint1(pos1)
            line.SetPoint2(pos2)
            line_mapper = vtk.vtkPolyDataMapper()
            line_mapper.SetInputConnection(line.GetOutputPort())
            line_actor = vtk.vtkActor()
            line_actor.SetMapper(line_mapper)
            self.renderer.AddActor(line_actor)

        self.renderer.SetBackground(0.1, 0.1, 0.1)
        self.renderer.ResetCamera()

    def snap_to_right(self):
        # Snap the window to the right side of the screen
        screen = QApplication.desktop().screenGeometry()
        screen_width = screen.width()
        screen_height = screen.height()

        window_width = self.width()
        window_height = self.height()

        self.setGeometry(QRect(screen_width - window_width, 0, window_width, window_height))


# Example usage in another script
if __name__ == "__main__":
    app = QApplication(sys.argv)
    viewer = Molecule3DViewer("CCO")  # Replace with your desired SMILES code
    sys.exit(app.exec_())
