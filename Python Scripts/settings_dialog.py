from PyQt5.QtWidgets import QDialog, QVBoxLayout, QCheckBox, QPushButton
from settings_manager import SettingsManager

class SettingsDialog(QDialog):
    def __init__(self, parent=None):
        super(SettingsDialog, self).__init__(parent)
        self.setWindowTitle("Settings")
        self.setGeometry(300, 300, 200, 200)
        self.settings_manager = SettingsManager()

        layout = QVBoxLayout()

        self.dark_mode_checkbox = QCheckBox("Enable Dark Mode", self)
        layout.addWidget(self.dark_mode_checkbox)

        self.dark_mode_checkbox.setChecked(parent.is_dark_mode)

        apply_button = QPushButton("Apply", self)
        apply_button.clicked.connect(self.apply_settings)
        layout.addWidget(apply_button)

        self.setLayout(layout)

    def apply_settings(self):
        self.parent().toggle_dark_mode(self.dark_mode_checkbox.isChecked())
        self.accept()
