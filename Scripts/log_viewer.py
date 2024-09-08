# log_viewer.py

import os
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QTextEdit, QFileDialog

class LogViewer(QDialog):
    """
    A class for viewing log files from the 'Logs' folder using a simple GUI.
    Automatically loads the 'molecule_viewer.log' file if available.
    """

    def __init__(self, logs_folder='Logs', log_file_name='molecule_viewer.log', parent=None):
        """
        Initialize the LogViewer, load the default log file ('molecule_viewer.log'), or allow the user to select one.

        Args:
            logs_folder (str): The folder where log files are stored (default is 'Logs').
            log_file_name (str): The default log file to load (default is 'molecule_viewer.log').
        """
        super().__init__(parent)

        # Set the root directory and log folder
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.logs_folder = os.path.join(root_dir, logs_folder)

        # Set up the GUI
        self.setWindowTitle("Log Viewer")
        self.setGeometry(300, 300, 600, 400)
        self.layout = QVBoxLayout()

        # Text area to display log content
        self.log_text_area = QTextEdit(self)
        self.log_text_area.setReadOnly(True)
        self.layout.addWidget(self.log_text_area)

        # Load the default log file if provided
        if log_file_name:
            self.load_log_file(log_file_name)

        # Buttons for selecting a log file and closing the viewer
        button_layout = QHBoxLayout()

        load_button = QPushButton("Load Log File", self)
        load_button.clicked.connect(self.select_log_file)
        button_layout.addWidget(load_button)

        close_button = QPushButton("Close", self)
        close_button.clicked.connect(self.close)
        button_layout.addWidget(close_button)

        self.layout.addLayout(button_layout)
        self.setLayout(self.layout)

    def load_log_file(self, file_name):
        """
        Load the content of a log file and display it in the text area.

        Args:
            file_name (str): The name of the log file to load.
        """
        log_file_path = os.path.join(self.logs_folder, file_name)

        if os.path.exists(log_file_path):
            with open(log_file_path, 'r') as file:
                self.log_text_area.setText(file.read())
        else:
            self.log_text_area.setText(f"Log file {file_name} not found.")

    def select_log_file(self):
        """
        Open a file dialog for the user to select a log file to view.
        """
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Log File", self.logs_folder, "Log Files (*.log *.txt)")
        if file_path:
            self.load_log_file(os.path.basename(file_path))
