import sys
from PyQt5.QtWidgets import QApplication
from viewer import MoleculeViewer

if __name__ == '__main__':
    app = QApplication(sys.argv)
    viewer = MoleculeViewer()
    viewer.show()
    sys.exit(app.exec_())

    #To Do:
    """
    Fix bond length calcs, add bond strength, add bond angles, 
    add 3D render, 
    add remove from database button, add sort functionality to database view, add dynamic columns configured in settings menu to display what user wants.
    Fix settings not saving, some functional groups are detected incorrectly (ketone detected with aldehyde, not detecting alkyne in C#C)
    

    """