from PySide6.QtWidgets import QApplication
import sys
from widgets import *


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
