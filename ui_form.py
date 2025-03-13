# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'form.ui'
##
## Created by: Qt User Interface Compiler version 6.8.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QAction, QBrush, QColor, QConicalGradient,
    QCursor, QFont, QFontDatabase, QGradient,
    QIcon, QImage, QKeySequence, QLinearGradient,
    QPainter, QPalette, QPixmap, QRadialGradient,
    QTransform)
from PySide6.QtWidgets import (QApplication, QButtonGroup, QCheckBox, QComboBox,
    QFrame, QGridLayout, QHBoxLayout, QLabel,
    QLineEdit, QMainWindow, QMenu, QMenuBar,
    QPushButton, QRadioButton, QSizePolicy, QSpacerItem,
    QStatusBar, QTabWidget, QVBoxLayout, QWidget)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.setEnabled(True)
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        icon = QIcon()
        icon.addFile(u"../alignment_calc_custom_pulse/icon.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        MainWindow.setWindowIcon(icon)
        self.saveTrace = QAction(MainWindow)
        self.saveTrace.setObjectName(u"saveTrace")
        self.saveTrace.setCheckable(False)
        self.saveTrace.setEnabled(False)
        self.actionExit = QAction(MainWindow)
        self.actionExit.setObjectName(u"actionExit")
        self.actionAbout = QAction(MainWindow)
        self.actionAbout.setObjectName(u"actionAbout")
        self.actionClose_figures = QAction(MainWindow)
        self.actionClose_figures.setObjectName(u"actionClose_figures")
        self.actionHelp = QAction(MainWindow)
        self.actionHelp.setObjectName(u"actionHelp")
        self.actionPrecalculate = QAction(MainWindow)
        self.actionPrecalculate.setObjectName(u"actionPrecalculate")
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        sizePolicy.setHeightForWidth(self.centralwidget.sizePolicy().hasHeightForWidth())
        self.centralwidget.setSizePolicy(sizePolicy)
        self.gridLayout = QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName(u"gridLayout")
        self.horizontalLayout_27 = QHBoxLayout()
        self.horizontalLayout_27.setObjectName(u"horizontalLayout_27")
        self.verticalLayout_9 = QVBoxLayout()
        self.verticalLayout_9.setObjectName(u"verticalLayout_9")
        self.label_4 = QLabel(self.centralwidget)
        self.label_4.setObjectName(u"label_4")

        self.verticalLayout_9.addWidget(self.label_4, 0, Qt.AlignmentFlag.AlignHCenter)

        self.InitConditionsTab = QTabWidget(self.centralwidget)
        self.InitConditionsTab.setObjectName(u"InitConditionsTab")
        self.InitConditionsTab.setTabShape(QTabWidget.TabShape.Rounded)
        self.Boltzman = QWidget()
        self.Boltzman.setObjectName(u"Boltzman")
        self.verticalLayout_5 = QVBoxLayout(self.Boltzman)
        self.verticalLayout_5.setObjectName(u"verticalLayout_5")
        self.horizontalLayout_16 = QHBoxLayout()
        self.horizontalLayout_16.setObjectName(u"horizontalLayout_16")
        self.label_6 = QLabel(self.Boltzman)
        self.label_6.setObjectName(u"label_6")
        self.label_6.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_16.addWidget(self.label_6)

        self.Temperature = QLineEdit(self.Boltzman)
        self.Temperature.setObjectName(u"Temperature")
        self.Temperature.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.Temperature.setMaxLength(8)
        self.Temperature.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_16.addWidget(self.Temperature)

        self.label_31 = QLabel(self.Boltzman)
        self.label_31.setObjectName(u"label_31")
        self.label_31.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_16.addWidget(self.label_31)


        self.verticalLayout_5.addLayout(self.horizontalLayout_16)

        self.horizontalLayout_17 = QHBoxLayout()
        self.horizontalLayout_17.setObjectName(u"horizontalLayout_17")
        self.label_32 = QLabel(self.Boltzman)
        self.label_32.setObjectName(u"label_32")
        self.label_32.setFrameShape(QFrame.Shape.NoFrame)
        self.label_32.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_17.addWidget(self.label_32)

        self.abundanceEven = QLineEdit(self.Boltzman)
        self.abundanceEven.setObjectName(u"abundanceEven")
        self.abundanceEven.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.abundanceEven.setMaxLength(8)
        self.abundanceEven.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_17.addWidget(self.abundanceEven)


        self.verticalLayout_5.addLayout(self.horizontalLayout_17)

        self.horizontalLayout_18 = QHBoxLayout()
        self.horizontalLayout_18.setObjectName(u"horizontalLayout_18")
        self.label_33 = QLabel(self.Boltzman)
        self.label_33.setObjectName(u"label_33")
        self.label_33.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_18.addWidget(self.label_33)

        self.abundanceOdd = QLineEdit(self.Boltzman)
        self.abundanceOdd.setObjectName(u"abundanceOdd")
        self.abundanceOdd.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.abundanceOdd.setMaxLength(8)
        self.abundanceOdd.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_18.addWidget(self.abundanceOdd)


        self.verticalLayout_5.addLayout(self.horizontalLayout_18)

        self.horizontalLayout_19 = QHBoxLayout()
        self.horizontalLayout_19.setObjectName(u"horizontalLayout_19")
        self.label_38 = QLabel(self.Boltzman)
        self.label_38.setObjectName(u"label_38")
        self.label_38.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_19.addWidget(self.label_38)

        self.num_ensemble_states = QLabel(self.Boltzman)
        self.num_ensemble_states.setObjectName(u"num_ensemble_states")
        self.num_ensemble_states.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_19.addWidget(self.num_ensemble_states)

        self.percentile = QLineEdit(self.Boltzman)
        self.percentile.setObjectName(u"percentile")
        self.percentile.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.percentile.setMaxLength(8)
        self.percentile.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_19.addWidget(self.percentile)

        self.label_40 = QLabel(self.Boltzman)
        self.label_40.setObjectName(u"label_40")
        self.label_40.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_19.addWidget(self.label_40)


        self.verticalLayout_5.addLayout(self.horizontalLayout_19)

        self.horizontalLayout_20 = QHBoxLayout()
        self.horizontalLayout_20.setObjectName(u"horizontalLayout_20")
        self.label_41 = QLabel(self.Boltzman)
        self.label_41.setObjectName(u"label_41")
        self.label_41.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_20.addWidget(self.label_41)

        self.ELfactor = QLineEdit(self.Boltzman)
        self.ELfactor.setObjectName(u"ELfactor")
        self.ELfactor.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.ELfactor.setMaxLength(8)
        self.ELfactor.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_20.addWidget(self.ELfactor)


        self.verticalLayout_5.addLayout(self.horizontalLayout_20)

        self.horizontalLayout_22 = QHBoxLayout()
        self.horizontalLayout_22.setObjectName(u"horizontalLayout_22")
        self.cos2d = QCheckBox(self.Boltzman)
        self.cos2d.setObjectName(u"cos2d")
        self.cos2d.setTristate(False)

        self.horizontalLayout_22.addWidget(self.cos2d)

        self.cos2dlabel = QLabel(self.Boltzman)
        self.cos2dlabel.setObjectName(u"cos2dlabel")
        self.cos2dlabel.setTextFormat(Qt.TextFormat.RichText)

        self.horizontalLayout_22.addWidget(self.cos2dlabel)

        self.horizontalSpacer_2 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_22.addItem(self.horizontalSpacer_2)


        self.verticalLayout_5.addLayout(self.horizontalLayout_22)

        self.label_5 = QLabel(self.Boltzman)
        self.label_5.setObjectName(u"label_5")
        sizePolicy.setHeightForWidth(self.label_5.sizePolicy().hasHeightForWidth())
        self.label_5.setSizePolicy(sizePolicy)

        self.verticalLayout_5.addWidget(self.label_5, 0, Qt.AlignmentFlag.AlignHCenter)

        self.InitConditionsTab.addTab(self.Boltzman, "")
        self.singleState = QWidget()
        self.singleState.setObjectName(u"singleState")
        self.verticalLayout_4 = QVBoxLayout(self.singleState)
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.label_34 = QLabel(self.singleState)
        self.label_34.setObjectName(u"label_34")
        font = QFont()
        font.setPointSize(18)
        font.setBold(True)
        self.label_34.setFont(font)

        self.verticalLayout_4.addWidget(self.label_34)

        self.horizontalLayout_21 = QHBoxLayout()
        self.horizontalLayout_21.setObjectName(u"horizontalLayout_21")
        self.label_35 = QLabel(self.singleState)
        self.label_35.setObjectName(u"label_35")
        self.label_35.setFont(font)

        self.horizontalLayout_21.addWidget(self.label_35)

        self.J = QLineEdit(self.singleState)
        self.J.setObjectName(u"J")
        self.J.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.J.setInputMethodHints(Qt.InputMethodHint.ImhNone)
        self.J.setMaxLength(8)
        self.J.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)
        self.J.setReadOnly(False)

        self.horizontalLayout_21.addWidget(self.J)

        self.K = QLineEdit(self.singleState)
        self.K.setObjectName(u"K")
        self.K.setEnabled(True)
        self.K.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.K.setMaxLength(8)
        self.K.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_21.addWidget(self.K)

        self.M = QLineEdit(self.singleState)
        self.M.setObjectName(u"M")
        self.M.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.M.setMaxLength(8)
        self.M.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_21.addWidget(self.M)

        self.label_36 = QLabel(self.singleState)
        self.label_36.setObjectName(u"label_36")
        self.label_36.setFont(font)

        self.horizontalLayout_21.addWidget(self.label_36)


        self.verticalLayout_4.addLayout(self.horizontalLayout_21)

        self.label_42 = QLabel(self.singleState)
        self.label_42.setObjectName(u"label_42")
        self.label_42.setTextFormat(Qt.TextFormat.AutoText)
        self.label_42.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft)

        self.verticalLayout_4.addWidget(self.label_42)

        self.InitConditionsTab.addTab(self.singleState, "")

        self.verticalLayout_9.addWidget(self.InitConditionsTab)


        self.horizontalLayout_27.addLayout(self.verticalLayout_9)


        self.gridLayout.addLayout(self.horizontalLayout_27, 1, 0, 1, 1)

        self.frame_3 = QFrame(self.centralwidget)
        self.frame_3.setObjectName(u"frame_3")
        self.frame_3.setFrameShape(QFrame.Shape.Panel)
        self.frame_3.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_3 = QVBoxLayout(self.frame_3)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.label_54 = QLabel(self.frame_3)
        self.label_54.setObjectName(u"label_54")

        self.verticalLayout_3.addWidget(self.label_54)

        self.horizontalLayout_12 = QHBoxLayout()
        self.horizontalLayout_12.setObjectName(u"horizontalLayout_12")
        self.label_39 = QLabel(self.frame_3)
        self.label_39.setObjectName(u"label_39")
        font1 = QFont()
        font1.setPointSize(12)
        font1.setBold(True)
        self.label_39.setFont(font1)

        self.horizontalLayout_12.addWidget(self.label_39)

        self.Jmax = QLineEdit(self.frame_3)
        self.Jmax.setObjectName(u"Jmax")
        self.Jmax.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.Jmax.setMaxLength(4)
        self.Jmax.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_12.addWidget(self.Jmax)


        self.verticalLayout_3.addLayout(self.horizontalLayout_12)

        self.label_55 = QLabel(self.frame_3)
        self.label_55.setObjectName(u"label_55")
        sizePolicy.setHeightForWidth(self.label_55.sizePolicy().hasHeightForWidth())
        self.label_55.setSizePolicy(sizePolicy)
        self.label_55.setMinimumSize(QSize(185, 31))
        self.label_55.setAlignment(Qt.AlignmentFlag.AlignLeading|Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignTop)
        self.label_55.setWordWrap(True)

        self.verticalLayout_3.addWidget(self.label_55)

        self.line_2 = QFrame(self.frame_3)
        self.line_2.setObjectName(u"line_2")
        self.line_2.setFrameShape(QFrame.Shape.HLine)
        self.line_2.setFrameShadow(QFrame.Shadow.Sunken)

        self.verticalLayout_3.addWidget(self.line_2)

        self.forceDT = QCheckBox(self.frame_3)
        self.forceDT.setObjectName(u"forceDT")

        self.verticalLayout_3.addWidget(self.forceDT)

        self.horizontalLayout_13 = QHBoxLayout()
        self.horizontalLayout_13.setObjectName(u"horizontalLayout_13")
        self.timestep = QLineEdit(self.frame_3)
        self.timestep.setObjectName(u"timestep")
        self.timestep.setEnabled(False)
        self.timestep.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_13.addWidget(self.timestep)

        self.label_28 = QLabel(self.frame_3)
        self.label_28.setObjectName(u"label_28")

        self.horizontalLayout_13.addWidget(self.label_28)


        self.verticalLayout_3.addLayout(self.horizontalLayout_13)

        self.choosePropTime = QCheckBox(self.frame_3)
        self.choosePropTime.setObjectName(u"choosePropTime")

        self.verticalLayout_3.addWidget(self.choosePropTime)

        self.horizontalLayout_14 = QHBoxLayout()
        self.horizontalLayout_14.setObjectName(u"horizontalLayout_14")
        self.propTime = QLineEdit(self.frame_3)
        self.propTime.setObjectName(u"propTime")
        self.propTime.setEnabled(False)
        self.propTime.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_14.addWidget(self.propTime)

        self.label_50 = QLabel(self.frame_3)
        self.label_50.setObjectName(u"label_50")

        self.horizontalLayout_14.addWidget(self.label_50)


        self.verticalLayout_3.addLayout(self.horizontalLayout_14)

        self.line = QFrame(self.frame_3)
        self.line.setObjectName(u"line")
        self.line.setFrameShape(QFrame.Shape.HLine)
        self.line.setFrameShadow(QFrame.Shadow.Sunken)

        self.verticalLayout_3.addWidget(self.line)

        self.label_45 = QLabel(self.frame_3)
        self.label_45.setObjectName(u"label_45")

        self.verticalLayout_3.addWidget(self.label_45)

        self.horizontalLayout_15 = QHBoxLayout()
        self.horizontalLayout_15.setObjectName(u"horizontalLayout_15")
        self.doODE = QRadioButton(self.frame_3)
        self.propagationGroup = QButtonGroup(MainWindow)
        self.propagationGroup.setObjectName(u"propagationGroup")
        self.propagationGroup.addButton(self.doODE)
        self.doODE.setObjectName(u"doODE")
        self.doODE.setChecked(True)

        self.horizontalLayout_15.addWidget(self.doODE)

        self.doMatrix = QRadioButton(self.frame_3)
        self.propagationGroup.addButton(self.doMatrix)
        self.doMatrix.setObjectName(u"doMatrix")

        self.horizontalLayout_15.addWidget(self.doMatrix)


        self.verticalLayout_3.addLayout(self.horizontalLayout_15)

        self.horizontalSpacer_3 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.verticalLayout_3.addItem(self.horizontalSpacer_3)

        self.mybutn = QPushButton(self.frame_3)
        self.mybutn.setObjectName(u"mybutn")
        self.mybutn.setCursor(QCursor(Qt.CursorShape.SizeVerCursor))

        self.verticalLayout_3.addWidget(self.mybutn)


        self.gridLayout.addWidget(self.frame_3, 1, 1, 1, 1)

        self.horizontalLayout_26 = QHBoxLayout()
        self.horizontalLayout_26.setSpacing(17)
        self.horizontalLayout_26.setObjectName(u"horizontalLayout_26")
        self.verticalLayout_8 = QVBoxLayout()
        self.verticalLayout_8.setObjectName(u"verticalLayout_8")
        self.frame = QFrame(self.centralwidget)
        self.frame.setObjectName(u"frame")
        self.frame.setFrameShape(QFrame.Shape.Panel)
        self.frame.setFrameShadow(QFrame.Shadow.Raised)
        self.frame.setLineWidth(1)
        self.verticalLayout = QVBoxLayout(self.frame)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.horizontalLayout_28 = QHBoxLayout()
        self.horizontalLayout_28.setObjectName(u"horizontalLayout_28")
        self.label = QLabel(self.frame)
        self.label.setObjectName(u"label")

        self.horizontalLayout_28.addWidget(self.label, 0, Qt.AlignmentFlag.AlignHCenter)


        self.verticalLayout.addLayout(self.horizontalLayout_28)

        self.label_11 = QLabel(self.frame)
        self.label_11.setObjectName(u"label_11")
        font2 = QFont()
        font2.setPointSize(12)
        self.label_11.setFont(font2)

        self.verticalLayout.addWidget(self.label_11)

        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.label_7 = QLabel(self.frame)
        self.label_7.setObjectName(u"label_7")
        self.label_7.setFont(font1)

        self.horizontalLayout.addWidget(self.label_7)

        self.Aconst = QLineEdit(self.frame)
        self.Aconst.setObjectName(u"Aconst")
        self.Aconst.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.Aconst.setMaxLength(11)
        self.Aconst.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout.addWidget(self.Aconst)

        self.label_9 = QLabel(self.frame)
        self.label_9.setObjectName(u"label_9")
        self.label_9.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing)

        self.horizontalLayout.addWidget(self.label_9)


        self.verticalLayout.addLayout(self.horizontalLayout)

        self.horizontalLayout_2 = QHBoxLayout()
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.label_8 = QLabel(self.frame)
        self.label_8.setObjectName(u"label_8")
        self.label_8.setFont(font1)

        self.horizontalLayout_2.addWidget(self.label_8)

        self.Bconst = QLineEdit(self.frame)
        self.Bconst.setObjectName(u"Bconst")
        self.Bconst.setMaxLength(11)
        self.Bconst.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_2.addWidget(self.Bconst)

        self.label_10 = QLabel(self.frame)
        self.label_10.setObjectName(u"label_10")
        self.label_10.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing)
        self.label_10.setWordWrap(False)

        self.horizontalLayout_2.addWidget(self.label_10)


        self.verticalLayout.addLayout(self.horizontalLayout_2)

        self.horizontalLayout_3 = QHBoxLayout()
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.label_51 = QLabel(self.frame)
        self.label_51.setObjectName(u"label_51")
        self.label_51.setFont(font1)

        self.horizontalLayout_3.addWidget(self.label_51)

        self.Dconst = QLineEdit(self.frame)
        self.Dconst.setObjectName(u"Dconst")
        self.Dconst.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.Dconst.setMaxLength(11)
        self.Dconst.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_3.addWidget(self.Dconst)

        self.label_53 = QLabel(self.frame)
        self.label_53.setObjectName(u"label_53")
        self.label_53.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing)

        self.horizontalLayout_3.addWidget(self.label_53)


        self.verticalLayout.addLayout(self.horizontalLayout_3)

        self.label_16 = QLabel(self.frame)
        self.label_16.setObjectName(u"label_16")
        self.label_16.setFont(font2)

        self.verticalLayout.addWidget(self.label_16)

        self.horizontalLayout_4 = QHBoxLayout()
        self.horizontalLayout_4.setObjectName(u"horizontalLayout_4")
        self.label_26 = QLabel(self.frame)
        self.label_26.setObjectName(u"label_26")
        self.label_26.setFont(font1)

        self.horizontalLayout_4.addWidget(self.label_26)

        self.alpha_par = QLineEdit(self.frame)
        self.alpha_par.setObjectName(u"alpha_par")
        self.alpha_par.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.alpha_par.setMaxLength(11)
        self.alpha_par.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_4.addWidget(self.alpha_par)

        self.label_23 = QLabel(self.frame)
        self.label_23.setObjectName(u"label_23")

        self.horizontalLayout_4.addWidget(self.label_23)


        self.verticalLayout.addLayout(self.horizontalLayout_4)

        self.horizontalLayout_5 = QHBoxLayout()
        self.horizontalLayout_5.setObjectName(u"horizontalLayout_5")
        self.label_25 = QLabel(self.frame)
        self.label_25.setObjectName(u"label_25")
        self.label_25.setFont(font1)

        self.horizontalLayout_5.addWidget(self.label_25)

        self.alpha_perp = QLineEdit(self.frame)
        self.alpha_perp.setObjectName(u"alpha_perp")
        self.alpha_perp.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.alpha_perp.setMaxLength(11)
        self.alpha_perp.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_5.addWidget(self.alpha_perp)

        self.label_24 = QLabel(self.frame)
        self.label_24.setObjectName(u"label_24")

        self.horizontalLayout_5.addWidget(self.label_24)


        self.verticalLayout.addLayout(self.horizontalLayout_5)

        self.horizontalLayout_6 = QHBoxLayout()
        self.horizontalLayout_6.setObjectName(u"horizontalLayout_6")
        self.label_27 = QLabel(self.frame)
        self.label_27.setObjectName(u"label_27")
        self.label_27.setFont(font1)

        self.horizontalLayout_6.addWidget(self.label_27)

        self.deltaAlpha = QLabel(self.frame)
        self.deltaAlpha.setObjectName(u"deltaAlpha")
        font3 = QFont()
        font3.setPointSize(11)
        font3.setBold(True)
        self.deltaAlpha.setFont(font3)
        self.deltaAlpha.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_6.addWidget(self.deltaAlpha)

        self.label_29 = QLabel(self.frame)
        self.label_29.setObjectName(u"label_29")

        self.horizontalLayout_6.addWidget(self.label_29)


        self.verticalLayout.addLayout(self.horizontalLayout_6)

        self.label_2 = QLabel(self.frame)
        self.label_2.setObjectName(u"label_2")

        self.verticalLayout.addWidget(self.label_2)

        self.moleculesBox = QComboBox(self.frame)
        self.moleculesBox.addItem("")
        self.moleculesBox.setObjectName(u"moleculesBox")
        self.moleculesBox.setInsertPolicy(QComboBox.InsertPolicy.InsertAlphabetically)

        self.verticalLayout.addWidget(self.moleculesBox)


        self.verticalLayout_8.addWidget(self.frame)


        self.horizontalLayout_26.addLayout(self.verticalLayout_8)


        self.gridLayout.addLayout(self.horizontalLayout_26, 0, 0, 1, 1)

        self.verticalLayout_7 = QVBoxLayout()
        self.verticalLayout_7.setObjectName(u"verticalLayout_7")
        self.PulseTab = QTabWidget(self.centralwidget)
        self.PulseTab.setObjectName(u"PulseTab")
        self.GaussianPulse = QWidget()
        self.GaussianPulse.setObjectName(u"GaussianPulse")
        self.verticalLayout_2 = QVBoxLayout(self.GaussianPulse)
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.label_3 = QLabel(self.GaussianPulse)
        self.label_3.setObjectName(u"label_3")

        self.verticalLayout_2.addWidget(self.label_3)

        self.horizontalLayout_7 = QHBoxLayout()
        self.horizontalLayout_7.setObjectName(u"horizontalLayout_7")
        self.label_22 = QLabel(self.GaussianPulse)
        self.label_22.setObjectName(u"label_22")

        self.horizontalLayout_7.addWidget(self.label_22)

        self.pulseEnergy = QLineEdit(self.GaussianPulse)
        self.pulseEnergy.setObjectName(u"pulseEnergy")
        self.pulseEnergy.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_7.addWidget(self.pulseEnergy)

        self.label_46 = QLabel(self.GaussianPulse)
        self.label_46.setObjectName(u"label_46")
        self.label_46.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing)

        self.horizontalLayout_7.addWidget(self.label_46)


        self.verticalLayout_2.addLayout(self.horizontalLayout_7)

        self.horizontalLayout_8 = QHBoxLayout()
        self.horizontalLayout_8.setObjectName(u"horizontalLayout_8")
        self.label_12 = QLabel(self.GaussianPulse)
        self.label_12.setObjectName(u"label_12")

        self.horizontalLayout_8.addWidget(self.label_12)

        self.pulseDuration = QLineEdit(self.GaussianPulse)
        self.pulseDuration.setObjectName(u"pulseDuration")
        self.pulseDuration.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_8.addWidget(self.pulseDuration)

        self.label_13 = QLabel(self.GaussianPulse)
        self.label_13.setObjectName(u"label_13")
        self.label_13.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing)

        self.horizontalLayout_8.addWidget(self.label_13)


        self.verticalLayout_2.addLayout(self.horizontalLayout_8)

        self.horizontalLayout_9 = QHBoxLayout()
        self.horizontalLayout_9.setObjectName(u"horizontalLayout_9")
        self.label_17 = QLabel(self.GaussianPulse)
        self.label_17.setObjectName(u"label_17")

        self.horizontalLayout_9.addWidget(self.label_17)

        self.pumpWaist = QLineEdit(self.GaussianPulse)
        self.pumpWaist.setObjectName(u"pumpWaist")
        self.pumpWaist.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_9.addWidget(self.pumpWaist)

        self.label_18 = QLabel(self.GaussianPulse)
        self.label_18.setObjectName(u"label_18")
        self.label_18.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing)

        self.horizontalLayout_9.addWidget(self.label_18)


        self.verticalLayout_2.addLayout(self.horizontalLayout_9)

        self.horizontalLayout_10 = QHBoxLayout()
        self.horizontalLayout_10.setObjectName(u"horizontalLayout_10")
        self.label_14 = QLabel(self.GaussianPulse)
        self.label_14.setObjectName(u"label_14")

        self.horizontalLayout_10.addWidget(self.label_14)

        self.pulseIntensity = QLineEdit(self.GaussianPulse)
        self.pulseIntensity.setObjectName(u"pulseIntensity")
        self.pulseIntensity.setEnabled(False)
        self.pulseIntensity.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_10.addWidget(self.pulseIntensity)

        self.label_15 = QLabel(self.GaussianPulse)
        self.label_15.setObjectName(u"label_15")
        self.label_15.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing)

        self.horizontalLayout_10.addWidget(self.label_15)


        self.verticalLayout_2.addLayout(self.horizontalLayout_10)

        self.horizontalLayout_11 = QHBoxLayout()
        self.horizontalLayout_11.setObjectName(u"horizontalLayout_11")
        self.label_47 = QLabel(self.GaussianPulse)
        self.label_47.setObjectName(u"label_47")

        self.horizontalLayout_11.addWidget(self.label_47)

        self.pulseFluence = QLineEdit(self.GaussianPulse)
        self.pulseFluence.setObjectName(u"pulseFluence")
        self.pulseFluence.setEnabled(False)
        self.pulseFluence.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_11.addWidget(self.pulseFluence)

        self.label_48 = QLabel(self.GaussianPulse)
        self.label_48.setObjectName(u"label_48")
        self.label_48.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing)

        self.horizontalLayout_11.addWidget(self.label_48)


        self.verticalLayout_2.addLayout(self.horizontalLayout_11)

        self.PulseTab.addTab(self.GaussianPulse, "")
        self.CustomPulse = QWidget()
        self.CustomPulse.setObjectName(u"CustomPulse")
        self.verticalLayout_10 = QVBoxLayout(self.CustomPulse)
        self.verticalLayout_10.setObjectName(u"verticalLayout_10")
        self.label_37 = QLabel(self.CustomPulse)
        self.label_37.setObjectName(u"label_37")

        self.verticalLayout_10.addWidget(self.label_37)

        self.horizontalLayout_29 = QHBoxLayout()
        self.horizontalLayout_29.setObjectName(u"horizontalLayout_29")
        self.label_52 = QLabel(self.CustomPulse)
        self.label_52.setObjectName(u"label_52")

        self.horizontalLayout_29.addWidget(self.label_52)

        self.pulseEnergy_2 = QLineEdit(self.CustomPulse)
        self.pulseEnergy_2.setObjectName(u"pulseEnergy_2")
        self.pulseEnergy_2.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_29.addWidget(self.pulseEnergy_2)

        self.label_61 = QLabel(self.CustomPulse)
        self.label_61.setObjectName(u"label_61")
        self.label_61.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing)

        self.horizontalLayout_29.addWidget(self.label_61)


        self.verticalLayout_10.addLayout(self.horizontalLayout_29)

        self.horizontalLayout_30 = QHBoxLayout()
        self.horizontalLayout_30.setObjectName(u"horizontalLayout_30")
        self.label_62 = QLabel(self.CustomPulse)
        self.label_62.setObjectName(u"label_62")

        self.horizontalLayout_30.addWidget(self.label_62)

        self.pumpWaist_2 = QLineEdit(self.CustomPulse)
        self.pumpWaist_2.setObjectName(u"pumpWaist_2")
        self.pumpWaist_2.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_30.addWidget(self.pumpWaist_2)

        self.label_59 = QLabel(self.CustomPulse)
        self.label_59.setObjectName(u"label_59")
        self.label_59.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing)

        self.horizontalLayout_30.addWidget(self.label_59)


        self.verticalLayout_10.addLayout(self.horizontalLayout_30)

        self.horizontalLayout_31 = QHBoxLayout()
        self.horizontalLayout_31.setObjectName(u"horizontalLayout_31")
        self.label_49 = QLabel(self.CustomPulse)
        self.label_49.setObjectName(u"label_49")

        self.horizontalLayout_31.addWidget(self.label_49)

        self.pulseFluence_2 = QLineEdit(self.CustomPulse)
        self.pulseFluence_2.setObjectName(u"pulseFluence_2")
        self.pulseFluence_2.setEnabled(False)
        self.pulseFluence_2.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_31.addWidget(self.pulseFluence_2)

        self.label_57 = QLabel(self.CustomPulse)
        self.label_57.setObjectName(u"label_57")
        self.label_57.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing)

        self.horizontalLayout_31.addWidget(self.label_57)


        self.verticalLayout_10.addLayout(self.horizontalLayout_31)

        self.loadpulse = QPushButton(self.CustomPulse)
        self.loadpulse.setObjectName(u"loadpulse")

        self.verticalLayout_10.addWidget(self.loadpulse)

        self.currentFile = QLineEdit(self.CustomPulse)
        self.currentFile.setObjectName(u"currentFile")
        self.currentFile.setEnabled(False)
        self.currentFile.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.verticalLayout_10.addWidget(self.currentFile)

        self.PulseTab.addTab(self.CustomPulse, "")

        self.verticalLayout_7.addWidget(self.PulseTab)

        self.frame_2 = QFrame(self.centralwidget)
        self.frame_2.setObjectName(u"frame_2")
        self.frame_2.setEnabled(True)
        sizePolicy.setHeightForWidth(self.frame_2.sizePolicy().hasHeightForWidth())
        self.frame_2.setSizePolicy(sizePolicy)
        self.frame_2.setFrameShape(QFrame.Shape.StyledPanel)
        self.frame_2.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_6 = QVBoxLayout(self.frame_2)
        self.verticalLayout_6.setObjectName(u"verticalLayout_6")
        self.horizontalLayout_23 = QHBoxLayout()
        self.horizontalLayout_23.setObjectName(u"horizontalLayout_23")
        self.label_43 = QLabel(self.frame_2)
        self.label_43.setObjectName(u"label_43")

        self.horizontalLayout_23.addWidget(self.label_43)

        self.t0 = QLineEdit(self.frame_2)
        self.t0.setObjectName(u"t0")
        self.t0.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_23.addWidget(self.t0)

        self.label_44 = QLabel(self.frame_2)
        self.label_44.setObjectName(u"label_44")
        self.label_44.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing)

        self.horizontalLayout_23.addWidget(self.label_44)


        self.verticalLayout_6.addLayout(self.horizontalLayout_23)

        self.horizontalLayout_24 = QHBoxLayout()
        self.horizontalLayout_24.setObjectName(u"horizontalLayout_24")
        self.label_19 = QLabel(self.frame_2)
        self.label_19.setObjectName(u"label_19")

        self.horizontalLayout_24.addWidget(self.label_19)

        self.probeWaist = QLineEdit(self.frame_2)
        self.probeWaist.setObjectName(u"probeWaist")
        self.probeWaist.setMaxLength(8)
        self.probeWaist.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_24.addWidget(self.probeWaist)

        self.label_20 = QLabel(self.frame_2)
        self.label_20.setObjectName(u"label_20")
        self.label_20.setAlignment(Qt.AlignmentFlag.AlignBottom|Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing)

        self.horizontalLayout_24.addWidget(self.label_20)


        self.verticalLayout_6.addLayout(self.horizontalLayout_24)

        self.label_21 = QLabel(self.frame_2)
        self.label_21.setObjectName(u"label_21")

        self.verticalLayout_6.addWidget(self.label_21)

        self.horizontalLayout_25 = QHBoxLayout()
        self.horizontalLayout_25.setObjectName(u"horizontalLayout_25")
        self.label_30 = QLabel(self.frame_2)
        self.label_30.setObjectName(u"label_30")

        self.horizontalLayout_25.addWidget(self.label_30)

        self.Nshells = QLineEdit(self.frame_2)
        self.Nshells.setObjectName(u"Nshells")
        self.Nshells.setMaxLength(3)
        self.Nshells.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignTrailing|Qt.AlignmentFlag.AlignVCenter)

        self.horizontalLayout_25.addWidget(self.Nshells)


        self.verticalLayout_6.addLayout(self.horizontalLayout_25)


        self.verticalLayout_7.addWidget(self.frame_2)


        self.gridLayout.addLayout(self.verticalLayout_7, 0, 1, 1, 1)

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 473, 21))
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        self.menuAbout = QMenu(self.menubar)
        self.menuAbout.setObjectName(u"menuAbout")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(u"statusbar")
        MainWindow.setStatusBar(self.statusbar)
        QWidget.setTabOrder(self.moleculesBox, self.Aconst)
        QWidget.setTabOrder(self.Aconst, self.Bconst)
        QWidget.setTabOrder(self.Bconst, self.alpha_par)
        QWidget.setTabOrder(self.alpha_par, self.alpha_perp)
        QWidget.setTabOrder(self.alpha_perp, self.pulseEnergy)
        QWidget.setTabOrder(self.pulseEnergy, self.pulseDuration)
        QWidget.setTabOrder(self.pulseDuration, self.pumpWaist)
        QWidget.setTabOrder(self.pumpWaist, self.t0)
        QWidget.setTabOrder(self.t0, self.pulseEnergy_2)
        QWidget.setTabOrder(self.pulseEnergy_2, self.pumpWaist_2)
        QWidget.setTabOrder(self.pumpWaist_2, self.probeWaist)
        QWidget.setTabOrder(self.probeWaist, self.Nshells)
        QWidget.setTabOrder(self.Nshells, self.Temperature)
        QWidget.setTabOrder(self.Temperature, self.abundanceEven)
        QWidget.setTabOrder(self.abundanceEven, self.abundanceOdd)
        QWidget.setTabOrder(self.abundanceOdd, self.percentile)
        QWidget.setTabOrder(self.percentile, self.ELfactor)
        QWidget.setTabOrder(self.ELfactor, self.J)
        QWidget.setTabOrder(self.J, self.K)
        QWidget.setTabOrder(self.K, self.M)
        QWidget.setTabOrder(self.M, self.Jmax)
        QWidget.setTabOrder(self.Jmax, self.cos2d)
        QWidget.setTabOrder(self.cos2d, self.forceDT)
        QWidget.setTabOrder(self.forceDT, self.timestep)
        QWidget.setTabOrder(self.timestep, self.doODE)
        QWidget.setTabOrder(self.doODE, self.doMatrix)
        QWidget.setTabOrder(self.doMatrix, self.mybutn)
        QWidget.setTabOrder(self.mybutn, self.PulseTab)
        QWidget.setTabOrder(self.PulseTab, self.InitConditionsTab)
        QWidget.setTabOrder(self.InitConditionsTab, self.pulseFluence)
        QWidget.setTabOrder(self.pulseFluence, self.pulseFluence_2)
        QWidget.setTabOrder(self.pulseFluence_2, self.loadpulse)
        QWidget.setTabOrder(self.loadpulse, self.currentFile)
        QWidget.setTabOrder(self.currentFile, self.pulseIntensity)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuAbout.menuAction())
        self.menuFile.addAction(self.saveTrace)
        self.menuFile.addAction(self.actionPrecalculate)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionExit)
        self.menuAbout.addAction(self.actionHelp)
        self.menuAbout.addAction(self.actionClose_figures)
        self.menuAbout.addSeparator()
        self.menuAbout.addAction(self.actionAbout)

        self.retranslateUi(MainWindow)

        self.InitConditionsTab.setCurrentIndex(0)
        self.PulseTab.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"Laser induced alignment calculator for symmetric top molecules", None))
        self.saveTrace.setText(QCoreApplication.translate("MainWindow", u"Save last trace data", None))
#if QT_CONFIG(shortcut)
        self.saveTrace.setShortcut(QCoreApplication.translate("MainWindow", u"Ctrl+S", None))
#endif // QT_CONFIG(shortcut)
        self.actionExit.setText(QCoreApplication.translate("MainWindow", u"Exit", None))
#if QT_CONFIG(shortcut)
        self.actionExit.setShortcut(QCoreApplication.translate("MainWindow", u"Esc", None))
#endif // QT_CONFIG(shortcut)
        self.actionAbout.setText(QCoreApplication.translate("MainWindow", u"About", None))
        self.actionClose_figures.setText(QCoreApplication.translate("MainWindow", u"Close figures", None))
#if QT_CONFIG(shortcut)
        self.actionClose_figures.setShortcut(QCoreApplication.translate("MainWindow", u"Q", None))
#endif // QT_CONFIG(shortcut)
        self.actionHelp.setText(QCoreApplication.translate("MainWindow", u"Help", None))
#if QT_CONFIG(shortcut)
        self.actionHelp.setShortcut(QCoreApplication.translate("MainWindow", u"?", None))
#endif // QT_CONFIG(shortcut)
        self.actionPrecalculate.setText(QCoreApplication.translate("MainWindow", u"Precalculate 2d matrix elements", None))
#if QT_CONFIG(shortcut)
        self.actionPrecalculate.setShortcut(QCoreApplication.translate("MainWindow", u"C", None))
#endif // QT_CONFIG(shortcut)
        self.label_4.setText(QCoreApplication.translate("MainWindow", u"Initial state", None))
        self.label_6.setText(QCoreApplication.translate("MainWindow", u"Temperature", None))
        self.Temperature.setText("")
        self.Temperature.setPlaceholderText(QCoreApplication.translate("MainWindow", u"0", None))
        self.label_31.setText(QCoreApplication.translate("MainWindow", u"Kelvin", None))
#if QT_CONFIG(tooltip)
        self.label_32.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Due to spin statistics and inversion/reflection symmetry, some rotational states either don't exists or are more abundant in nature than others.<br/><br/>For example, the symmetric molecule CS2 can only occupy even J states, as <span style=\" vertical-align:super;\">32</span>S is a boson with spin 0 and the electronic ground state <span style=\" vertical-align:super;\">1</span><span style=\" font-weight:600;\">\u03a3</span><span style=\" font-weight:600; vertical-align:sub;\">g</span> is even. The total wavefunction should be even. It is electronic*rotational*nuclear. The nuclear and electronic parts are even, so the rotational wavefunction must also be even. O<span style=\" vertical-align:sub;\">2</span> can only exist in odd J states. Species with nonzero nuclear spin (e.g. H<span style=\" vertical-align:sub;\">2</span>) can have &quot;strange&quot; abundancies due to the difference between the number of odd and even nuclear states. Hydrogen for example is three times more likely to be "
                        "in an odd J state than an even.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_32.setText(QCoreApplication.translate("MainWindow", u"Even abundance", None))
        self.abundanceEven.setText("")
        self.abundanceEven.setPlaceholderText(QCoreApplication.translate("MainWindow", u"1", None))
#if QT_CONFIG(tooltip)
        self.label_33.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Due to spin statistics and inversion/reflection symmetry, some rotational states either don't exists or are more abundant in nature than others.<br/><br/>For example, the symmetric molecule CS2 can only occupy even J states, as <span style=\" vertical-align:super;\">32</span>S is a boson with spin 0 and the electronic ground state <span style=\" vertical-align:super;\">1</span><span style=\" font-weight:600;\">\u03a3</span><span style=\" font-weight:600; vertical-align:sub;\">g</span> is even. The total wavefunction should be even. It is electronic*rotational*nuclear. The nuclear and electronic parts are even, so the rotational wavefunction must also be even. O<span style=\" vertical-align:sub;\">2</span> can only exist in odd J states. Species with nonzero nuclear spin (e.g. H<span style=\" vertical-align:sub;\">2</span>) can have &quot;strange&quot; abundancies due to the difference between the number of odd and even nuclear states. Hydrogen for example is three times more likely to be "
                        "in an odd J state than an even.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_33.setText(QCoreApplication.translate("MainWindow", u"Odd abundance", None))
        self.abundanceOdd.setText("")
        self.abundanceOdd.setPlaceholderText(QCoreApplication.translate("MainWindow", u"1", None))
#if QT_CONFIG(tooltip)
        self.label_38.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Number of ensemble states to propagate to resolve 99.9% of the partition function.</p><p>Note: this can be much lower than the actual number of states in a thermal ensemble. However, due to symmetry, not all ensemble memebers needs to be propagated.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_38.setText(QCoreApplication.translate("MainWindow", u"Ensemble size", None))
#if QT_CONFIG(tooltip)
        self.num_ensemble_states.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Multiply with #FVA shells to get number of initial states to solve the Schr\u00f6dinger equation for.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.num_ensemble_states.setText(QCoreApplication.translate("MainWindow", u"N/A", None))
#if QT_CONFIG(tooltip)
        self.percentile.setToolTip(QCoreApplication.translate("MainWindow", u"Resolve this percentile of the Boltzmann distribution", None))
#endif // QT_CONFIG(tooltip)
        self.percentile.setText("")
        self.percentile.setPlaceholderText(QCoreApplication.translate("MainWindow", u"99.9", None))
#if QT_CONFIG(tooltip)
        self.label_40.setToolTip(QCoreApplication.translate("MainWindow", u"Resolve this percentile of the Boltzmann distribution", None))
#endif // QT_CONFIG(tooltip)
        self.label_40.setText(QCoreApplication.translate("MainWindow", u"%", None))
#if QT_CONFIG(tooltip)
        self.label_41.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>The Even Lavie valve tends to supress low M states, such that the equi-partion theorem is violated. To model this, weigh each state by this factor raised to the |M|'th power. If the factor is smaller than 1, it suppresses high M states. If it is larger than 1, low M states are suppressed.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_41.setText(QCoreApplication.translate("MainWindow", u"Even Lavie M anisotropy", None))
#if QT_CONFIG(tooltip)
        self.ELfactor.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>The Even Lavie valve tends to supress low M states, such that the equi-partion theorem is violated. To model this, weigh each state by this factor raised to the |M|'th power. If the factor is smaller than 1, it suppresses high M states. If it is larger than 1, low M states are suppressed.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.ELfactor.setText("")
        self.ELfactor.setPlaceholderText(QCoreApplication.translate("MainWindow", u"1.0", None))
#if QT_CONFIG(tooltip)
        self.cos2d.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>&lt;cos<span style=\" vertical-align:super;\">2</span>\u03b8<span style=\" vertical-align:sub;\">2d</span>&gt; is typically the only available observable in an actual experiment.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.cos2d.setStatusTip(QCoreApplication.translate("MainWindow", u"For direct comparison with experiments", None))
#endif // QT_CONFIG(statustip)
        self.cos2d.setText("")
#if QT_CONFIG(tooltip)
        self.cos2dlabel.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>&lt;cos<span style=\" vertical-align:super;\">2</span>\u03b8<span style=\" vertical-align:sub;\">2d</span>&gt; is typically the only available observable in an actual experiment.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.cos2dlabel.setStatusTip(QCoreApplication.translate("MainWindow", u"For direct comparison with experiments", None))
#endif // QT_CONFIG(statustip)
        self.cos2dlabel.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Also calculate &lt;cos<span style=\" vertical-align:super;\">2</span> \u03b8<span style=\" vertical-align:sub;\">2d</span>&gt;</p></body></html>", None))
#if QT_CONFIG(tooltip)
        self.label_5.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>If you don't, the matrix elements will be re-computed every time you calculate a trace.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_5.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Tip: In the file menu, you can calculate<br/>&lt;J'KM|cos<span style=\" vertical-align:super;\">2</span>\u03b8<span style=\" vertical-align:sub;\">2D</span>|JKM&gt; matrix elements in<br/>advance for faster calculations.</p></body></html>", None))
        self.InitConditionsTab.setTabText(self.InitConditionsTab.indexOf(self.Boltzman), QCoreApplication.translate("MainWindow", u"Boltzmann ensemble", None))
        self.label_34.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>|JKM&gt; = </p></body></html>", None))
        self.label_35.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>|</p></body></html>", None))
#if QT_CONFIG(tooltip)
        self.J.setToolTip(QCoreApplication.translate("MainWindow", u"J quantum number", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.J.setStatusTip(QCoreApplication.translate("MainWindow", u"Total angular momentum quantum number", None))
#endif // QT_CONFIG(statustip)
        self.J.setText("")
        self.J.setPlaceholderText(QCoreApplication.translate("MainWindow", u"0", None))
#if QT_CONFIG(tooltip)
        self.K.setToolTip(QCoreApplication.translate("MainWindow", u"K quantum number (always 0 for linear molecules)", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.K.setStatusTip(QCoreApplication.translate("MainWindow", u"Projection on molecular axis", None))
#endif // QT_CONFIG(statustip)
        self.K.setText("")
        self.K.setPlaceholderText(QCoreApplication.translate("MainWindow", u"0", None))
#if QT_CONFIG(tooltip)
        self.M.setToolTip(QCoreApplication.translate("MainWindow", u"M quantum number", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.M.setStatusTip(QCoreApplication.translate("MainWindow", u"Projection on laser axis", None))
#endif // QT_CONFIG(statustip)
        self.M.setText("")
        self.M.setPlaceholderText(QCoreApplication.translate("MainWindow", u"0", None))
        self.label_36.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>&gt;</p></body></html>", None))
#if QT_CONFIG(tooltip)
        self.label_42.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Without focal volume averaging, also indicate the distribution of populated J states with &lt;J&gt; and &lt;J&gt;\u00b1sqrt(&lt;J<span style=\" vertical-align:super;\">2</span>&gt;-&lt;J&gt;<span style=\" vertical-align:super;\">2</span>).</p><p>The J value J<span style=\" vertical-align:sub;\">99.9%</span>, below which 99.9% of the probability amplitude is found, is also plotted.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_42.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>If #FVA shells is 1, also plot extra<br/>information about the wave function.<br/>(hover mouse here for more information)</p></body></html>", None))
        self.InitConditionsTab.setTabText(self.InitConditionsTab.indexOf(self.singleState), QCoreApplication.translate("MainWindow", u"Single |JKM> state", None))
        self.label_54.setText(QCoreApplication.translate("MainWindow", u"Basis size, propagation (test)", None))
#if QT_CONFIG(tooltip)
        self.label_39.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>The lower, the faster the calculation. If too low, the result will be incorrect.</p><p>The optimal value depends on the laser pulse and the molecule in question.</p><p>Tip: for adiabatic alignment, Jmax can be really low.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_39.setText(QCoreApplication.translate("MainWindow", u"Jmax =", None))
#if QT_CONFIG(tooltip)
        self.Jmax.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>The lower, the faster the calculation. If too low, the result will be incorrect.</p><p>The optimal value depends on the laser pulse and the molecule in question.</p><p>Tip: for adiabatic alignment, Jmax can be really low.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.Jmax.setText("")
        self.Jmax.setPlaceholderText(QCoreApplication.translate("MainWindow", u"140", None))
#if QT_CONFIG(tooltip)
        self.label_55.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>The lower, the faster the calculation. If too low, the result will be incorrect.</p><p>The optimal value depends on the laser pulse and the molecule in question.</p><p>Tip: for adiabatic alignment, Jmax can be really low.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_55.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Largest J quantum number |JKM&gt;. Should be as low as possible.</p></body></html>", None))
#if QT_CONFIG(tooltip)
        self.forceDT.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Sampling time step, not propagation time step. By default, this corresponds to 5 times the Nyquist sampling rate for traces in the impulsive limit.</p><p>If your simulate adiabatic alignment, consider increasing the time step a lot to save time.</p><p>Note that for custom pulses, this determines the sampling</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.forceDT.setStatusTip(QCoreApplication.translate("MainWindow", u"For plotting only", None))
#endif // QT_CONFIG(statustip)
        self.forceDT.setText(QCoreApplication.translate("MainWindow", u"Force time step", None))
#if QT_CONFIG(tooltip)
        self.timestep.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Sampling time step, not propagation time step. By default, this corresponds to 5 times the Nyquist sampling rate for traces in the impulsive limit.</p><p>If your simulate adiabatic alignment, consider increasing the time step a lot to save time.</p><p>Note that for custom pulses, this determines the sampling</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_28.setText(QCoreApplication.translate("MainWindow", u"ps", None))
#if QT_CONFIG(tooltip)
        self.choosePropTime.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Time to propagate pulse to.</p><p>By default propagates to first revival</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.choosePropTime.setStatusTip("")
#endif // QT_CONFIG(statustip)
        self.choosePropTime.setText(QCoreApplication.translate("MainWindow", u"Choose propagation time", None))
#if QT_CONFIG(tooltip)
        self.propTime.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Time to propagate pulse to.</p><p>By default propagates to first revival</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_50.setText(QCoreApplication.translate("MainWindow", u"ps", None))
        self.label_45.setText(QCoreApplication.translate("MainWindow", u"Propagation method:", None))
#if QT_CONFIG(tooltip)
        self.doODE.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Solve Schr\u00f6dinger equation for the expansion coefficients of the wave function in the |JKM&gt; basis as an ordinary differential equation.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.doODE.setText(QCoreApplication.translate("MainWindow", u"ODE", None))
#if QT_CONFIG(tooltip)
        self.doMatrix.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Solve the Schr\u00f6dinger equation by applying the time evolution operator to the initial states. This is typically faster for very long pulses.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.doMatrix.setText(QCoreApplication.translate("MainWindow", u"Matrix", None))
#if QT_CONFIG(tooltip)
        self.mybutn.setToolTip(QCoreApplication.translate("MainWindow", u"Perform the calculation", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.mybutn.setStatusTip(QCoreApplication.translate("MainWindow", u"Click to calculate an alignment trace", None))
#endif // QT_CONFIG(statustip)
        self.mybutn.setText(QCoreApplication.translate("MainWindow", u"Go!", None))
        self.label.setText(QCoreApplication.translate("MainWindow", u"Molecule data", None))
        self.label_11.setText(QCoreApplication.translate("MainWindow", u"Rotational constants", None))
#if QT_CONFIG(tooltip)
        self.label_7.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Secondary rotational constant.</p><p>Linear molecules have A=0.<br/>Prolate tops have A&gt;B.<br/>Oblate tops have A&lt;B.<br/>Spherical tops have A=B</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.label_7.setStatusTip(QCoreApplication.translate("MainWindow", u"Linear molecules have A=0", None))
#endif // QT_CONFIG(statustip)
        self.label_7.setText(QCoreApplication.translate("MainWindow", u"A", None))
#if QT_CONFIG(tooltip)
        self.Aconst.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Secondary rotational constant.</p><p>Linear molecules have A=0.<br/>Prolate tops have A&gt;B.<br/>Oblate tops have A&lt;B.<br/>Spherical tops have A=B</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.Aconst.setStatusTip(QCoreApplication.translate("MainWindow", u"Linear molecules have A=0", None))
#endif // QT_CONFIG(statustip)
        self.Aconst.setText("")
        self.Aconst.setPlaceholderText(QCoreApplication.translate("MainWindow", u"0", None))
#if QT_CONFIG(tooltip)
        self.label_9.setToolTip(QCoreApplication.translate("MainWindow", u"Linear molecules have A=0", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.label_9.setStatusTip(QCoreApplication.translate("MainWindow", u"Linear molecules have A=0", None))
#endif // QT_CONFIG(statustip)
        self.label_9.setText(QCoreApplication.translate("MainWindow", u"GHz", None))
#if QT_CONFIG(tooltip)
        self.label_8.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Main rotational constant.</p><p>This rotational constant corresponds to the two identical moments of inertia in a principal frame.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_8.setText(QCoreApplication.translate("MainWindow", u"B", None))
#if QT_CONFIG(tooltip)
        self.Bconst.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Main rotational constant.</p><p>This rotational constant corresponds to the two identical moments of inertia in a principal frame.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.Bconst.setPlaceholderText("")
        self.label_10.setText(QCoreApplication.translate("MainWindow", u"GHz", None))
#if QT_CONFIG(tooltip)
        self.label_51.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Rotational distortion constant. Calculation is faster if this is zero</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.label_51.setStatusTip(QCoreApplication.translate("MainWindow", u"Linear molecules have A=0", None))
#endif // QT_CONFIG(statustip)
        self.label_51.setText(QCoreApplication.translate("MainWindow", u"D<sub>J</sub>", None))
#if QT_CONFIG(tooltip)
        self.Dconst.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Rotational distortion constant. Calculation is faster if this is zero</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.Dconst.setStatusTip(QCoreApplication.translate("MainWindow", u"Linear molecules have A=0", None))
#endif // QT_CONFIG(statustip)
        self.Dconst.setText("")
        self.Dconst.setPlaceholderText(QCoreApplication.translate("MainWindow", u"0", None))
#if QT_CONFIG(tooltip)
        self.label_53.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Rotational distortion constant. Calculation is faster if this is zero</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.label_53.setStatusTip(QCoreApplication.translate("MainWindow", u"Linear molecules have A=0", None))
#endif // QT_CONFIG(statustip)
        self.label_53.setText(QCoreApplication.translate("MainWindow", u"KHz", None))
        self.label_16.setText(QCoreApplication.translate("MainWindow", u"Polarizability tensor", None))
#if QT_CONFIG(tooltip)
        self.label_26.setToolTip(QCoreApplication.translate("MainWindow", u"Polarizability volume along main molecular axis", None))
#endif // QT_CONFIG(tooltip)
        self.label_26.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" font-size:12pt;\">\u03b1</span><span style=\" font-size:12pt; vertical-align:sub;\">\u2225</span></p></body></html>", None))
#if QT_CONFIG(tooltip)
        self.alpha_par.setToolTip(QCoreApplication.translate("MainWindow", u"Polarizability volume along main molecular axis", None))
#endif // QT_CONFIG(tooltip)
        self.alpha_par.setText("")
#if QT_CONFIG(tooltip)
        self.label_23.setToolTip("")
#endif // QT_CONFIG(tooltip)
        self.label_23.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>\u00c5<span style=\" vertical-align:super;\">3</span></p></body></html>", None))
#if QT_CONFIG(tooltip)
        self.label_25.setToolTip(QCoreApplication.translate("MainWindow", u"Polarizability volume perpendicular to the main molecular axis", None))
#endif // QT_CONFIG(tooltip)
        self.label_25.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>\u03b1<span style=\" vertical-align:sub;\">&perp;</span></p></body></html>", None))
#if QT_CONFIG(tooltip)
        self.alpha_perp.setToolTip(QCoreApplication.translate("MainWindow", u"Polarizability volume perpendicular to the main molecular axis", None))
#endif // QT_CONFIG(tooltip)
        self.alpha_perp.setText("")
        self.label_24.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>\u00c5<span style=\" vertical-align:super;\">3</span></p></body></html>", None))
#if QT_CONFIG(tooltip)
        self.label_27.setToolTip(QCoreApplication.translate("MainWindow", u"Polarizability volume anisotropy", None))
#endif // QT_CONFIG(tooltip)
        self.label_27.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>\u0394\u03b1</p></body></html>", None))
#if QT_CONFIG(tooltip)
        self.deltaAlpha.setToolTip(QCoreApplication.translate("MainWindow", u"Polarizability volume anisotropy", None))
#endif // QT_CONFIG(tooltip)
        self.deltaAlpha.setText(QCoreApplication.translate("MainWindow", u"0", None))
        self.label_29.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>\u00c5<span style=\" vertical-align:super;\">3</span></p></body></html>", None))
        self.label_2.setText(QCoreApplication.translate("MainWindow", u"Predefined molecules", None))
        self.moleculesBox.setItemText(0, QCoreApplication.translate("MainWindow", u"Load...", None))

#if QT_CONFIG(tooltip)
        self.moleculesBox.setToolTip("")
#endif // QT_CONFIG(tooltip)
        self.label_3.setText(QCoreApplication.translate("MainWindow", u"Gaussian, linearly polarized", None))
#if QT_CONFIG(tooltip)
        self.label_22.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Pulse duration, Full width at half max.</p><p>You are allowed to use scientific notation. E.g. setting the duration to 1e6 fs would set the pulse duration to 1 ns, and you would be simulating adiabatic alignment.</p><p>When doing so, it is very advantageous to lower the basis size drastically.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_22.setText(QCoreApplication.translate("MainWindow", u"Pulse Energy", None))
#if QT_CONFIG(tooltip)
        self.pulseEnergy.setToolTip(QCoreApplication.translate("MainWindow", u"Pump pulse energy", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.pulseEnergy.setStatusTip(QCoreApplication.translate("MainWindow", u"Pump pulse energy", None))
#endif // QT_CONFIG(statustip)
        self.pulseEnergy.setText("")
        self.pulseEnergy.setPlaceholderText(QCoreApplication.translate("MainWindow", u"100", None))
        self.label_46.setText(QCoreApplication.translate("MainWindow", u"\u00b5J", None))
#if QT_CONFIG(tooltip)
        self.label_12.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Pulse duration, Full width at half max.</p><p>You are allowed to use scientific notation. E.g. setting the duration to 1e6 fs would set the pulse duration to 1 ns, and you would be simulating adiabatic alignment.</p><p>When doing so, it is very advantageous to lower the basis size drastically.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_12.setText(QCoreApplication.translate("MainWindow", u"FWHM duration", None))
#if QT_CONFIG(tooltip)
        self.pulseDuration.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Pulse duration, Full width at half max.</p><p>You are allowed to use scientific notation. E.g. setting the duration to 1e6 fs would set the pulse duration to 1 ns, and you would be simulating adiabatic alignment.</p><p>When doing so, it is very advantageous to lower the basis size drastically.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.pulseDuration.setText("")
        self.pulseDuration.setPlaceholderText(QCoreApplication.translate("MainWindow", u"300", None))
        self.label_13.setText(QCoreApplication.translate("MainWindow", u"fs", None))
#if QT_CONFIG(tooltip)
        self.label_17.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Spot size of the pump (or kick) pulse</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_17.setText(QCoreApplication.translate("MainWindow", u"Pump waist", None))
#if QT_CONFIG(tooltip)
        self.pumpWaist.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Spot size of the pump (or kick) pulse</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.pumpWaist.setText("")
        self.pumpWaist.setPlaceholderText(QCoreApplication.translate("MainWindow", u"35", None))
        self.label_18.setText(QCoreApplication.translate("MainWindow", u"\u00b5m", None))
#if QT_CONFIG(tooltip)
        self.label_14.setToolTip(QCoreApplication.translate("MainWindow", u"Kick (or pump) pulse", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.label_14.setStatusTip(QCoreApplication.translate("MainWindow", u"Kick pulse intensity", None))
#endif // QT_CONFIG(statustip)
        self.label_14.setText(QCoreApplication.translate("MainWindow", u"Peak intensity", None))
#if QT_CONFIG(tooltip)
        self.pulseIntensity.setToolTip(QCoreApplication.translate("MainWindow", u"Kick (or pump) pulse", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.pulseIntensity.setStatusTip(QCoreApplication.translate("MainWindow", u"Kick pulse intensity", None))
#endif // QT_CONFIG(statustip)
        self.pulseIntensity.setText("")
        self.pulseIntensity.setPlaceholderText(QCoreApplication.translate("MainWindow", u"16.3", None))
        self.label_15.setText(QCoreApplication.translate("MainWindow", u"TW/cm\u00b2", None))
#if QT_CONFIG(tooltip)
        self.label_47.setToolTip(QCoreApplication.translate("MainWindow", u"Kick (or pump) pulse", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.label_47.setStatusTip(QCoreApplication.translate("MainWindow", u"Kick pulse intensity", None))
#endif // QT_CONFIG(statustip)
        self.label_47.setText(QCoreApplication.translate("MainWindow", u"Fluence", None))
#if QT_CONFIG(tooltip)
        self.pulseFluence.setToolTip(QCoreApplication.translate("MainWindow", u"Kick (or pump) pulse", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.pulseFluence.setStatusTip(QCoreApplication.translate("MainWindow", u"Kick pulse intensity", None))
#endif // QT_CONFIG(statustip)
        self.pulseFluence.setText("")
        self.pulseFluence.setPlaceholderText(QCoreApplication.translate("MainWindow", u"5.2", None))
        self.label_48.setText(QCoreApplication.translate("MainWindow", u"J/cm\u00b2", None))
        self.PulseTab.setTabText(self.PulseTab.indexOf(self.GaussianPulse), QCoreApplication.translate("MainWindow", u"Gaussian Pulse", None))
        self.label_37.setText(QCoreApplication.translate("MainWindow", u"Custom Pulse, linearly polarized", None))
#if QT_CONFIG(tooltip)
        self.label_52.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Pulse duration, Full width at half max.</p><p>You are allowed to use scientific notation. E.g. setting the duration to 1e6 fs would set the pulse duration to 1 ns, and you would be simulating adiabatic alignment.</p><p>When doing so, it is very advantageous to lower the basis size drastically.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_52.setText(QCoreApplication.translate("MainWindow", u"Pulse Energy", None))
#if QT_CONFIG(tooltip)
        self.pulseEnergy_2.setToolTip(QCoreApplication.translate("MainWindow", u"Kick (or pump) pulse", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.pulseEnergy_2.setStatusTip(QCoreApplication.translate("MainWindow", u"Kick pulse intensity", None))
#endif // QT_CONFIG(statustip)
        self.pulseEnergy_2.setText("")
        self.pulseEnergy_2.setPlaceholderText(QCoreApplication.translate("MainWindow", u"100", None))
        self.label_61.setText(QCoreApplication.translate("MainWindow", u"\u00b5J", None))
#if QT_CONFIG(tooltip)
        self.label_62.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Spot size of the pump (or kick) pulse</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_62.setText(QCoreApplication.translate("MainWindow", u"Pump waist", None))
#if QT_CONFIG(tooltip)
        self.pumpWaist_2.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Spot size of the pump (or kick) pulse</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.pumpWaist_2.setText("")
        self.pumpWaist_2.setPlaceholderText(QCoreApplication.translate("MainWindow", u"35", None))
        self.label_59.setText(QCoreApplication.translate("MainWindow", u"\u00b5m", None))
#if QT_CONFIG(tooltip)
        self.label_49.setToolTip(QCoreApplication.translate("MainWindow", u"Kick (or pump) pulse", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.label_49.setStatusTip(QCoreApplication.translate("MainWindow", u"Kick pulse intensity", None))
#endif // QT_CONFIG(statustip)
        self.label_49.setText(QCoreApplication.translate("MainWindow", u"Fluence", None))
#if QT_CONFIG(tooltip)
        self.pulseFluence_2.setToolTip(QCoreApplication.translate("MainWindow", u"Kick (or pump) pulse", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.pulseFluence_2.setStatusTip(QCoreApplication.translate("MainWindow", u"Kick pulse intensity", None))
#endif // QT_CONFIG(statustip)
        self.pulseFluence_2.setText("")
        self.pulseFluence_2.setPlaceholderText(QCoreApplication.translate("MainWindow", u"5.2", None))
        self.label_57.setText(QCoreApplication.translate("MainWindow", u"J/cm\u00b2", None))
        self.loadpulse.setText(QCoreApplication.translate("MainWindow", u"Load Pulse", None))
#if QT_CONFIG(tooltip)
        self.currentFile.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Spot size of the pump (or kick) pulse</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.currentFile.setText("")
        self.currentFile.setPlaceholderText("")
        self.PulseTab.setTabText(self.PulseTab.indexOf(self.CustomPulse), QCoreApplication.translate("MainWindow", u"Custom Pulse", None))
#if QT_CONFIG(tooltip)
        self.label_43.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Kick pulse time</p><p>Tip: If you need to propagate to longer times, add a second, 0 intensity pulse at t = t0 + the extra time you need.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_43.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>T<span style=\" vertical-align:sub;\">0</span></p></body></html>", None))
#if QT_CONFIG(tooltip)
        self.t0.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Kick pulse time</p><p>Tip: If you need to propagate to longer times, add a second, 0 intensity pulse at t = t0 + the extra time you need.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.t0.setText("")
        self.t0.setPlaceholderText(QCoreApplication.translate("MainWindow", u"0", None))
        self.label_44.setText(QCoreApplication.translate("MainWindow", u"ps", None))
#if QT_CONFIG(tooltip)
        self.label_19.setToolTip(QCoreApplication.translate("MainWindow", u"Spot size of the probe pulse", None))
#endif // QT_CONFIG(tooltip)
        self.label_19.setText(QCoreApplication.translate("MainWindow", u"Probe waist", None))
#if QT_CONFIG(tooltip)
        self.probeWaist.setToolTip(QCoreApplication.translate("MainWindow", u"Spot size of the probe pulse", None))
#endif // QT_CONFIG(tooltip)
        self.probeWaist.setText("")
        self.probeWaist.setPlaceholderText(QCoreApplication.translate("MainWindow", u"25", None))
        self.label_20.setText(QCoreApplication.translate("MainWindow", u"\u00b5m", None))
        self.label_21.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" text-decoration: underline;\">F</span>ocal <span style=\" text-decoration: underline;\">V</span>olume <span style=\" text-decoration: underline;\">A</span>veraging shell(s)</p></body></html>", None))
#if QT_CONFIG(tooltip)
        self.label_30.setToolTip(QCoreApplication.translate("MainWindow", u"Number of focal volume iso intensity shells to average traces over", None))
#endif // QT_CONFIG(tooltip)
        self.label_30.setText(QCoreApplication.translate("MainWindow", u"#FVA shells", None))
#if QT_CONFIG(tooltip)
        self.Nshells.setToolTip(QCoreApplication.translate("MainWindow", u"Number of focal volume iso intensity shells to average traces over", None))
#endif // QT_CONFIG(tooltip)
        self.Nshells.setText("")
        self.Nshells.setPlaceholderText(QCoreApplication.translate("MainWindow", u"1", None))
        self.menuFile.setTitle(QCoreApplication.translate("MainWindow", u"File", None))
        self.menuAbout.setTitle(QCoreApplication.translate("MainWindow", u"Help", None))
    # retranslateUi

