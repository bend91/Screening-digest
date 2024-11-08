from PySide6.QtWidgets import QApplication
import sys
from widgets import *


# TODO - Add a list of band sizes for each sequence, maybe so that can be copied?


class MainApplication(QApplication):
    def __init__(self, *args, **kwargs):
        QApplication.__init__(self)


def main():
    main_app = MainApplication(sys.argv)
    window = MainWidget()
    window.show()
    main_app.exec()

if __name__ == "__main__":
    main()
