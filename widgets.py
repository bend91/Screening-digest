from PySide6.QtWidgets import QWidget, QTextEdit, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QRadioButton, QScrollArea, QGraphicsTextItem
from PySide6.QtCharts import QChart, QScatterSeries, QValueAxis, QLogValueAxis, QChartView
from PySide6.QtGui import QPainterPath, QPainter, QImage, QBrush, QPen
from PySide6.QtCore import QRectF, QSize, Qt, QPointF
from functions import *


# TODO - make this class only accept dna bases
class DNATextEdit(QTextEdit):
    def __init__(self, *args, **kwargs):
        QTextEdit.__init__(self, *args, **kwargs)


class GelMarkerSeries(QScatterSeries):
    def __init__(self, *args, **kwargs):
        QScatterSeries.__init__(self, *args, **kwargs)
        marker_path = QPainterPath()
        marker_path.addRect(QRectF(0, 0, 20, 3))
        image = QImage(QSize(100, 100), QImage.Format_ARGB32)
        painter = QPainter(image)
        # painter.fillRect(image.rect(), Qt.white)
        painter.fillPath(marker_path, Qt.black)
        painter.end()
        brush = QBrush(image)
        self.setMarkerSize(20)
        self.setMarkerShape(QScatterSeries.MarkerShapeRectangle)
        self.setBrush(brush)
        self.setPen(QPen(Qt.transparent))
        self.setPointLabelsVisible(True)
        self.setPointLabelsColor(Qt.black)
        self.setPointLabelsFormat("@yPoint")


class LadderGelMarkerSeries(GelMarkerSeries):
    def __init__(self, *args, **kwargs):
        GelMarkerSeries.__init__(self)
        self.setBrush(QBrush(Qt.transparent))
        self.setPen(QPen(Qt.transparent))
        self.setPointLabelsVisible(True)
        self.setPointLabelsColor(Qt.black)
        self.setPointLabelsFormat("@yPoint")


class GelPlotChart(QChart):
    def __init__(self, *args, **kwargs):
        QChart.__init__(self, *args, **kwargs)
        self.x_axis = QValueAxis()
        self.y_axis = QLogValueAxis()

        self.x_axis.setGridLineVisible(False)
        self.y_axis.setGridLineVisible(False)

        self.x_axis.setRange(-1, 3)
        self.y_axis.setRange(100, 15000)
        self.setAxisX(self.x_axis)
        self.setAxisY(self.y_axis)


class MainWidget(QWidget):
    def __init__(self, *args, **kwargs):
        QWidget.__init__(self)

        self.enzymes_dict = {}
        self.rbutton_dict = {}
        self.selected_enzyme = None


        layout = QVBoxLayout()
        h_layout = QHBoxLayout()
        v_layout1 = QVBoxLayout()
        v_layout2 = QVBoxLayout()
        label1 = QLabel("Enter Sequence 1")
        label2 = QLabel("Enter Sequence 2")

        self.text_box1 = DNATextEdit()
        self.text_box2 = DNATextEdit()
        self.results_layout = QHBoxLayout()
        self.results_chart = GelPlotChart()
        self.results_view = QChartView(self.results_chart)
        self.results_scroll_area = QScrollArea()
        self.results_layout.addWidget(self.results_scroll_area)
        self.results_layout.addWidget(self.results_view)

        v_layout1.addWidget(label1)
        v_layout1.addWidget(self.text_box1)

        v_layout2.addWidget(label2)
        v_layout2.addWidget(self.text_box2)

        self.submit_button = QPushButton("Submit")
        self.submit_button.clicked.connect(self.submit_clicked)

        h_layout.addLayout(v_layout1)
        h_layout.addLayout(v_layout2)

        layout.addLayout(h_layout)
        layout.addWidget(self.submit_button)
        layout.addLayout(self.results_layout)
        self.setLayout(layout)
        self.setGeometry(500, 500, 1000, 1000)

    def submit_clicked(self, *args, **kwargs):
        sequence1 = self.text_box1.toPlainText()
        sequence2 = self.text_box2.toPlainText()
        c1, c2, nc = sequence_comparison(sequence1, sequence2)
        best_cutters = identify_best_cutter(sequence1, sequence2, nc)

        enzyme_list_widget = QWidget()
        enz_layout = QVBoxLayout()

        # test = QTextEdit()
        # enz_layout.addWidget(test)
        # enz_layout.addWidget(QPushButton("test"))
        # enz_layout.addWidget(RadioButton("test"))

        for enzyme in best_cutters:
            print(enzyme)
            self.rbutton_dict[enzyme[0]] = RadioButton(enzyme[0])
            self.rbutton_dict[enzyme[0]].clicked.connect(self.radio_button_clicked)
            self.enzymes_dict[enzyme[0]] = (False, enzyme[1], enzyme[2])
            enz_layout.addWidget(self.rbutton_dict[enzyme[0]])
        # enzyme_list_widget.setLayout(enz_layout)
        self.results_layout.addLayout(enz_layout)
        # self.results_scroll_area.setWidget(enzyme_list_widget)


    def radio_button_clicked(self, *args, **kwargs):
        for k, v in self.rbutton_dict.items():
            if v.isChecked():
                self.selected_enzyme = v.text()
        self.plot_gel()

    def plot_gel(self, *args, **kwargs):
        ladder_choice = "HyperLadderI"
        ladder_plotted = False
        for series in self.results_chart.series():
            if ladder_choice in series.name():
                ladder_plotted = True
            else:
                self.results_chart.removeSeries(series)
        if not ladder_plotted:
            ladder = list(ladders.loc[ladder_choice].values)
            series0 = LadderGelMarkerSeries()
            series0.setName(f"{ladder_choice}_labels")
            series1 = GelMarkerSeries()
            series1.setName(ladder_choice)
            for i in ladder:
                series0.append(-0.5, i)
                series1.append(0, i)
            for series in [series0, series1]:
                self.results_chart.addSeries(series)
                series.attachAxis(self.results_chart.x_axis)
                series.attachAxis(self.results_chart.y_axis)

        series2 = GelMarkerSeries()
        series2.setName("Sequence 1")
        series3 = GelMarkerSeries()
        series3.setName("Sequence 2")

        data_series = [series2, series3]

        for i, series in enumerate(data_series, 1):
            series.append([QPointF(i, x) for x in self.enzymes_dict[self.selected_enzyme][i]])

        for series in data_series:
            self.results_chart.addSeries(series)
            series.attachAxis(self.results_chart.x_axis)
            series.attachAxis(self.results_chart.y_axis)
        ...

# /Users/benjamindraper/Library/CloudStorage/Dropbox/01_UCL/Archive/Plasmids/[PL00279]pCCL[GGStuff].SFFV.RFP.dna
# /Users/benjamindraper/Library/CloudStorage/Dropbox/01_UCL/Archive/Plasmids/[PL00303]pCCL.TE9_BBz-2A-19B11[L-H].dna
# /Users/benjamindraper/Library/CloudStorage/Dropbox/01_UCL/Archive/Plasmids/[PL00304]pCCL.TE9_BBz-2A-33G7[H-L].dna
# /Users/benjamindraper/Library/CloudStorage/Dropbox/01_UCL/Archive/Plasmids/[PL00305]pCCL.TE9_BBz-2A-33G7[L-H].dna
# /Users/benjamindraper/Library/CloudStorage/Dropbox/01_UCL/Archive/Plasmids/[PL00306]pCCL.TE9_BBz-2A-NbE10.dna
# /Users/benjamindraper/Library/CloudStorage/Dropbox/01_UCL/Archive/Plasmids/[PL00307]pCCL.TE9_BBz-2A-NbE10_ABD094.dna
# /Users/benjamindraper/Library/CloudStorage/Dropbox/01_UCL/Archive/Plasmids/[PL00308]pCCL.TE9_BBz-2A-sCD74.dna
# /Users/benjamindraper/Library/CloudStorage/Dropbox/01_UCL/Archive/Plasmids/[PL00225]pLKO.MIFsgRNA3.H2B-2A-mCherry.dna
# NotI, PstI, PstI

class RadioButton(QRadioButton):
    def __init__(self, *args, **kwargs):
        QRadioButton.__init__(self, *args, **kwargs)
        self.name = self.text()
        self.clicked.connect(self.selected)

    def selected(self, *args, **kwargs):
        # return self.text()
        print(self.text())

