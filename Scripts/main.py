import sys
from PyQt5.QtWidgets import QApplication
from viewer import MoleculeViewer

if __name__ == '__main__':
    app = QApplication(sys.argv)
    viewer = MoleculeViewer()
    viewer.show()
    sys.exit(app.exec_())
