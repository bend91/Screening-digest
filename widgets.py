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

        for enzyme in best_cutters:
            self.rbutton_dict[enzyme[0]] = RadioButton(enzyme[0])
            self.rbutton_dict[enzyme[0]].clicked.connect(self.radio_button_clicked)
            self.enzymes_dict[enzyme[0]] = (False, enzyme[1], enzyme[2])
            enz_layout.addWidget(self.rbutton_dict[enzyme[0]])
        enzyme_list_widget.setLayout(enz_layout)
        scroll_area = QScrollArea()
        scroll_area.setWidget(enzyme_list_widget)
        self.results_layout.addWidget(scroll_area)
        self.results_layout.addWidget(self.results_view)
        # self.qplot_gel()

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


class RadioButton(QRadioButton):
    def __init__(self, *args, **kwargs):
        QRadioButton.__init__(self, *args, **kwargs)
        self.name = self.text()
        self.clicked.connect(self.selected)

    def selected(self, *args, **kwargs):
        # return self.text()
        print(self.text())

### TextEdit options
"""
a.acceptDrops(                  a.focusInEvent(                 a.lower(                        a.setCornerWidget(              a.setWindowModified(
a.acceptRichText(               a.focusNextChild(               a.mapFrom(                      a.setCurrentCharFormat(         a.setWindowOpacity(
a.accessibleDescription(        a.focusNextPrevChild(           a.mapFromGlobal(                a.setCurrentFont(               a.setWindowRole(
a.accessibleName(               a.focusOutEvent(                a.mapFromParent(                a.setCursor(                    a.setWindowState(
a.actionEvent(                  a.focusPolicy(                  a.mapTo(                        a.setCursorWidth(               a.setWindowTitle(
a.actions(                      a.focusPreviousChild(           a.mapToGlobal(                  a.setDisabled(                  a.setWordWrapMode(
a.activateWindow(               a.focusProxy(                   a.mapToParent(                  a.setDocument(                  a.Shadow(
a.addAction(                    a.focusWidget(                  a.mask(                         a.setDocumentTitle(             a.Shape(
a.addActions(                   a.font(                         a.maximumHeight(                a.setEnabled(                   a.sharedPainter(
a.addScrollBarWidget(           a.fontFamily(                   a.maximumSize(                  a.setExtraSelections(           a.show(
a.adjustSize(                   a.fontInfo(                     a.maximumViewportSize(          a.setFixedHeight(               a.showEvent(
a.alignment(                    a.fontItalic(                   a.maximumWidth(                 a.setFixedSize(                 a.showFullScreen(
a.anchorAt(                     a.fontMetrics(                  a.mergeCurrentCharFormat(       a.setFixedWidth(                a.showMaximized(
a.append(                       a.fontPointSize(                a.metaObject(                   a.setFocus(                     a.showMinimized(
a.autoFillBackground(           a.fontUnderline(                a.metric(                       a.setFocusPolicy(               a.showNormal(
a.autoFormatting(               a.fontWeight(                   a.midLineWidth(                 a.setFocusProxy(                a.signalsBlocked(
a.AutoFormattingFlag(           a.foregroundRole(               a.minimumHeight(                a.setFont(                      a.size(
a.backgroundRole(               a.frameGeometry(                a.minimumSize(                  a.setFontFamily(                a.sizeAdjustPolicy(
a.backingStore(                 a.frameRect(                    a.minimumSizeHint(              a.setFontItalic(                a.SizeAdjustPolicy(
a.baseSize(                     a.frameShadow(                  a.minimumWidth(                 a.setFontPointSize(             a.sizeHint(
a.blockSignals(                 a.frameShape(                   a.mouseDoubleClickEvent(        a.setFontUnderline(             a.sizeIncrement(
a.canInsertFromMimeData(        a.frameSize(                    a.mouseGrabber(                 a.setFontWeight(                a.sizePolicy(
a.canPaste(                     a.frameStyle(                   a.mouseMoveEvent(               a.setForegroundRole(            a.stackUnder(
a.changeEvent(                  a.frameWidth(                   a.mousePressEvent(              a.setFrameRect(                 a.startTimer(
a.childAt(                      a.geometry(                     a.mouseReleaseEvent(            a.setFrameShadow(               a.staticMetaObject
a.childEvent(                   a.grab(                         a.move(                         a.setFrameShape(                a.statusTip(
a.children(                     a.grabGesture(                  a.moveCursor(                   a.setFrameStyle(                a.style(
a.childrenRect(                 a.grabKeyboard(                 a.moveEvent(                    a.setGeometry(                  a.StyleMask(
a.childrenRegion(               a.grabMouse(                    a.moveToThread(                 a.setGraphicsEffect(            a.styleSheet(
a.clear(                        a.grabShortcut(                 a.nativeEvent(                  a.setHidden(                    a.tabChangesFocus(
a.clearFocus(                   a.graphicsEffect(               a.nativeParentWidget(           a.setHorizontalScrollBar(       a.tabletEvent(
a.clearMask(                    a.graphicsProxyWidget(          a.nextInFocusChain(             a.setHorizontalScrollBarPolicy( a.tabStopDistance(
a.close(                        a.hasFocus(                     a.normalGeometry(               a.setHtml(                      a.testAttribute(
a.closeEvent(                   a.hasHeightForWidth(            a.objectName(                   a.setInputMethodHints(          a.textBackgroundColor(
a.colorCount(                   a.hasMouseTracking(             a.objectNameChanged(            a.setLayout(                    a.textChanged(
a.connect(                      a.hasTabletTracking(            a.overrideWindowFlags(          a.setLayoutDirection(           a.textColor(
a.connectNotify(                a.height(                       a.overrideWindowState(          a.setLineWidth(                 a.textCursor(
a.contentsMargins(              a.heightForWidth(               a.overwriteMode(                a.setLineWrapColumnOrWidth(     a.textInteractionFlags(
a.contentsRect(                 a.heightMM(                     a.PaintDeviceMetric(            a.setLineWrapMode(              a.thread(
a.contextMenuEvent(             a.hide(                         a.paintEngine(                  a.setLocale(                    a.timerEvent(
a.contextMenuPolicy(            a.hideEvent(                    a.painters                      a.setMarkdown(                  a.toHtml(
a.copy(                         a.horizontalScrollBar(          a.paintEvent(                   a.setMask(                      a.toMarkdown(
a.copyAvailable(                a.horizontalScrollBarPolicy(    a.paintingActive(               a.setMaximumHeight(             a.toolTip(
a.cornerWidget(                 a.inherits(                     a.palette(                      a.setMaximumSize(               a.toolTipDuration(
a.create(                       a.initPainter(                  a.parent(                       a.setMaximumWidth(              a.toPlainText(
a.createMimeDataFromSelection(  a.initStyleOption(              a.parentWidget(                 a.setMidLineWidth(              a.topLevelWidget(
a.createStandardContextMenu(    a.inputMethodEvent(             a.paste(                        a.setMinimumHeight(             a.tr(
a.createWindowContainer(        a.inputMethodHints(             a.physicalDpiX(                 a.setMinimumSize(               a.underMouse(
a.createWinId(                  a.inputMethodQuery(             a.physicalDpiY(                 a.setMinimumWidth(              a.undo(
a.currentCharFormat(            a.insertAction(                 a.placeholderText(              a.setMouseTracking(             a.undoAvailable(
a.currentCharFormatChanged(     a.insertActions(                a.pos(                          a.setObjectName(                a.ungrabGesture(
a.currentFont(                  a.insertFromMimeData(           a.previousInFocusChain(         a.setOverwriteMode(             a.unsetCursor(
a.cursor(                       a.insertHtml(                   a.print_(                       a.setPalette(                   a.unsetLayoutDirection(
a.cursorForPosition(            a.insertPlainText(              a.property(                     a.setParent(                    a.unsetLocale(
a.cursorPositionChanged(        a.installEventFilter(           a.raise_(                       a.setPlaceholderText(           a.update(
a.cursorRect(                   a.internalWinId(                a.receivers(                    a.setPlainText(                 a.updateGeometry(
a.cursorWidth(                  a.isActiveWindow(               a.rect(                         a.setProperty(                  a.updateMicroFocus(
a.customContextMenuRequested(   a.isAncestorOf(                 a.redirected(                   a.setReadOnly(                  a.updatesEnabled(
a.customEvent(                  a.isEnabled(                    a.redo(                         a.setScreen(                    a.verticalScrollBar(
a.cut(                          a.isEnabledTo(                  a.redoAvailable(                a.setShortcutAutoRepeat(        a.verticalScrollBarPolicy(
a.deleteLater(                  a.isFullScreen(                 a.releaseKeyboard(              a.setShortcutEnabled(           a.viewport(
a.depth(                        a.isHidden(                     a.releaseMouse(                 a.setSizeAdjustPolicy(          a.viewportEvent(
a.destroy(                      a.isLeftToRight(                a.releaseShortcut(              a.setSizeIncrement(             a.viewportMargins(
a.destroyed(                    a.isMaximized(                  a.removeAction(                 a.setSizePolicy(                a.viewportSizeHint(
a.devicePixelRatio(             a.isMinimized(                  a.removeEventFilter(            a.setStatusTip(                 a.visibleRegion(
a.devicePixelRatioF(            a.isModal(                      a.render(                       a.setStyle(                     a.whatsThis(
a.devicePixelRatioFScale(       a.isQuickItemType(              a.RenderFlag(                   a.setStyleSheet(                a.wheelEvent(
a.devType(                      a.isReadOnly(                   a.repaint(                      a.setTabChangesFocus(           a.width(
a.disconnect(                   a.isRightToLeft(                a.resize(                       a.setTabletTracking(            a.widthMM(
a.disconnectNotify(             a.isSignalConnected(            a.resizeEvent(                  a.setTabOrder(                  a.window(
a.document(                     a.isTopLevel(                   a.restoreGeometry(              a.setTabStopDistance(           a.windowFilePath(
a.documentTitle(                a.isUndoRedoEnabled(            a.saveGeometry(                 a.setText(                      a.windowFlags(
a.doSetTextCursor(              a.isVisible(                    a.screen(                       a.setTextBackgroundColor(       a.windowHandle(
a.dragEnterEvent(               a.isVisibleTo(                  a.scroll(                       a.setTextColor(                 a.windowIcon(
a.dragLeaveEvent(               a.isWidgetType(                 a.scrollBarWidgets(             a.setTextCursor(                a.windowIconChanged(
a.dragMoveEvent(                a.isWindow(                     a.scrollContentsBy(             a.setTextInteractionFlags(      a.windowIconText(
a.drawFrame(                    a.isWindowModified(             a.scrollToAnchor(               a.setToolTip(                   a.windowIconTextChanged(
a.dropEvent(                    a.isWindowType(                 a.selectAll(                    a.setToolTipDuration(           a.windowModality(
a.dumpObjectInfo(               a.keyboardGrabber(              a.selectionChanged(             a.setUndoRedoEnabled(           a.windowOpacity(
a.dumpObjectTree(               a.keyPressEvent(                a.sender(                       a.setUpdatesEnabled(            a.windowRole(
a.dynamicPropertyNames(         a.keyReleaseEvent(              a.senderSignalIndex(            a.setupViewport(                a.windowState(
a.effectiveWinId(               a.killTimer(                    a.setAcceptDrops(               a.setVerticalScrollBar(         a.windowTitle(
a.emit(                         a.layout(                       a.setAcceptRichText(            a.setVerticalScrollBarPolicy(   a.windowTitleChanged(
a.ensureCursorVisible(          a.layoutDirection(              a.setAccessibleDescription(     a.setViewport(                  a.windowType(
a.ensurePolished(               a.leaveEvent(                   a.setAccessibleName(            a.setViewportMargins(           a.winId(
a.enterEvent(                   a.lineWidth(                    a.setAlignment(                 a.setVisible(                   a.wordWrapMode(
a.event(                        a.lineWrapColumnOrWidth(        a.setAttribute(                 a.setWhatsThis(                 a.x(
a.eventFilter(                  a.lineWrapMode(                 a.setAutoFillBackground(        a.setWindowFilePath(            a.y(
a.ExtraSelection(               a.LineWrapMode(                 a.setAutoFormatting(            a.setWindowFlag(                a.zoomIn(
a.extraSelections(              a.loadResource(                 a.setBackgroundRole(            a.setWindowFlags(               a.zoomInF(
a.find(                         a.locale(                       a.setBaseSize(                  a.setWindowIcon(                a.zoomOut(
a.findChild(                    a.logicalDpiX(                  a.setContentsMargins(           a.setWindowIconText(
a.findChildren(                 a.logicalDpiY(                  a.setContextMenuPolicy(         a.setWindowModality(
"""
