import os
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QTextEdit

class LogViewer(QDialog):
    """
    A generic log viewer class that displays the contents of a log file in a dialog.
    
    Attributes:
    log_folder (str): The folder where the log file is located.
    log_file_name (str): The name of the log file to display.
    """
    
    def __init__(self, log_folder, log_file_name):
        """
        Initialize the LogViewer dialog.
        
        Args:
        log_folder (str): Path to the folder where the log file is stored.
        log_file_name (str): Name of the log file to display.
        """
        super().__init__()
        self.log_folder = log_folder
        self.log_file_name = log_file_name
        self.setWindowTitle("Log Viewer")  # Set the title of the dialog
        self.setGeometry(100, 100, 600, 400)  # Set the size and position of the window

        # Create a layout to hold the log text display
        layout = QVBoxLayout()

        # Create a QTextEdit widget for displaying the log file content (read-only)
        self.log_text_edit = QTextEdit(self)
        self.log_text_edit.setReadOnly(True)  # Make the text box read-only
        layout.addWidget(self.log_text_edit)  # Add the text box to the layout

        self.setLayout(layout)  # Set the layout for the dialog
        self.load_log_file()  # Load the log file content on initialization

    def load_log_file(self):
        """
        Load the log file and display its content in the QTextEdit widget.
        """
        log_file_path = os.path.join(self.log_folder, self.log_file_name)  # Construct full log file path
        
        try:
            # Try to open and read the log file content
            with open(log_file_path, 'r') as log_file:
                self.log_text_edit.setPlainText(log_file.read())  # Display the content in the text box
        except FileNotFoundError:
            # If the log file is not found, display an error message in the text box
            self.log_text_edit.setPlainText("Log file not found.")
