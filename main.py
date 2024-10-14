from PySide6.QtWidgets import QApplication, QWidget, QTextEdit, QVBoxLayout, QHBoxLayout, QPushButton
import sys
sys.path.append("/Users/benjamindraper/GitHub/ULA/Plasmid_Editor")
from ScreeningDigest_tk.functions import *

# Only needed for access to command line arguments
import sys

# You need one (and only one) QApplication instance per application.
# Pass in sys.argv to allow command line arguments for your app.
# If you know you won't use command line arguments QApplication([]) works too.


def MainApplication(QApplication):
    def __init__(self, *args, **kwargs):
        super(self, QApplication).__init__(*args, **kwargs)


def MainWindow(QWidget):
    def __init__(self, *args, **kwargs):
        # QWidget.__init__(self, *args, **kwargs)
        layout = QVBoxLayout()
        h_layout = QHBoxLayout()

        text_box1 = QTextEdit()
        text_box2 = QTextEdit()
        for widget in [text_box1, text_box2]:
            h_layout.addWidget(widget)
        submit_button = QPushButton("Submit")
        layout.addLayout(h_layout)
        layout.addWidget(submit_button)
        self.setLayout(layout)


def main():
    main_app = MainApplication(sys.argv)

# Create a Qt widget, which will be our window.
    window = MainWindow()
    # layout = QVBoxLayout()
    # h_layout = QHBoxLayout()
    
    # text_box1 = QTextEdit()
    # text_box2 = QTextEdit()

    # h_layout.addWidget(text_box1)
    # h_layout.addWidget(text_box2)

    # submit_button = QPushButton("Submit")

    # layout.addWidget(text_box)
    # window.setLayout(layout)
    window.show()
    main_app.exec()


if __name__ == "__main__":
    main()

