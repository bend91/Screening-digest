from PySide6.QtWidgets import QApplication

# Only needed for access to command line arguments
import sys
from widgets import *

# You need one (and only one) QApplication instance per application.
# Pass in sys.argv to allow command line arguments for your app.
# If you know you won't use command line arguments QApplication([]) works too.


class MainApplication(QApplication):
    def __init__(self, *args, **kwargs):
        QApplication.__init__(self)
        # super(self, QApplication).__init__(*args, **kwargs)


def main():
    main_app = MainApplication(sys.argv)
    # main_app = QApplication(sys.argv)


    window = MainWidget()
    # window = QWidget()
    # layout = QVBoxLayout()
    # h_layout = QHBoxLayout()
    # v_layout1 = QVBoxLayout()
    # v_layout2 = QVBoxLayout()
    
    # label1 = QLabel("Enter Sequence 1")
    # label2 = QLabel("Enter Sequence 2")
    # text_box1 = QTextEdit()
    # text_box2 = QTextEdit()

    # v_layout1.addWidget(label1)
    # v_layout1.addWidget(text_box1)

    # v_layout2.addWidget(label2)
    # v_layout2.addWidget(text_box2)

    # submit_button = QPushButton("Submit")

    # h_layout.addLayout(v_layout1)
    # h_layout.addLayout(v_layout2)

    # layout.addLayout(h_layout)
    # layout.addWidget(submit_button)
    # window.setLayout(layout)
    window.show()
    main_app.exec()

if __name__ == "__main__":
    main()

