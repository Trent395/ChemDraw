from PyQt5 import QtGui, QtCore

class DarkMode:
    @staticmethod
    def apply(app):
        dark_palette = QtGui.QPalette()

        # General colors for dark mode
        dark_palette.setColor(QtGui.QPalette.Window, QtGui.QColor(53, 53, 53))
        dark_palette.setColor(QtGui.QPalette.WindowText, QtCore.Qt.white)
        dark_palette.setColor(QtGui.QPalette.Base, QtGui.QColor(25, 25, 25))
        dark_palette.setColor(QtGui.QPalette.AlternateBase, QtGui.QColor(53, 53, 53))
        dark_palette.setColor(QtGui.QPalette.ToolTipBase, QtCore.Qt.white)
        dark_palette.setColor(QtGui.QPalette.ToolTipText, QtCore.Qt.black)
        dark_palette.setColor(QtGui.QPalette.Text, QtCore.Qt.white)

        # Button colors
        dark_palette.setColor(QtGui.QPalette.Button, QtGui.QColor(75, 75, 75))
        dark_palette.setColor(QtGui.QPalette.ButtonText, QtCore.Qt.white)

        # Input field colors
        dark_palette.setColor(QtGui.QPalette.Text, QtCore.Qt.white)
        dark_palette.setColor(QtGui.QPalette.PlaceholderText, QtGui.QColor(169, 169, 169))

        # Highlight colors
        dark_palette.setColor(QtGui.QPalette.Highlight, QtGui.QColor(42, 130, 218))
        dark_palette.setColor(QtGui.QPalette.HighlightedText, QtCore.Qt.white)

        app.setPalette(dark_palette)
