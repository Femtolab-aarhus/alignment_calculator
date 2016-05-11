# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'aboutForm.ui'
#
# Created: Wed May 11 09:14:47 2016
#      by: PyQt5 UI code generator 5.3.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(504, 254)
        self.label = QtWidgets.QLabel(Form)
        self.label.setGeometry(QtCore.QRect(12, 10, 120, 90))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(Form)
        self.label_2.setGeometry(QtCore.QRect(160, 10, 331, 201))
        self.label_2.setTextFormat(QtCore.Qt.RichText)
        self.label_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.label_2.setOpenExternalLinks(True)
        self.label_2.setObjectName("label_2")
        self.label_5 = QtWidgets.QLabel(Form)
        self.label_5.setGeometry(QtCore.QRect(10, 110, 141, 16))
        self.label_5.setTextFormat(QtCore.Qt.RichText)
        self.label_5.setOpenExternalLinks(True)
        self.label_5.setObjectName("label_5")
        self.Ok = QtWidgets.QPushButton(Form)
        self.Ok.setGeometry(QtCore.QRect(210, 220, 91, 24))
        self.Ok.setObjectName("Ok")

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "About"))
        self.label.setText(_translate("Form", "<html><head/><body><img src=\"./smalllogo.png\"/></body></html>"))
        self.label_2.setText(_translate("Form", "<html><head/><body><p>Laser induced alignment trace calculator.</p><p>Copyright 2016 Anders A. Søndergaard,<br/>PhD student in <a href=\"http://femtolab.au.dk\"><span style=\" text-decoration: underline; color:#0057ae;\">Femtolab</span></a> 2013-2016.<br/>Supervisor: Henrik Stapelfeldt.<br/><br/>This program calculates alignment traces for<br/>symmetric top molecules exposed to ultra-short<br/>laser pulses.</p><p>This program comes with ABSOLUTELY NO WARRANTY.<br/>This is free software, and you are welcome to<br/>redistribute it under certain conditions;<br/>see the LICENSE file for details.</p></body></html>"))
        self.label_5.setText(_translate("Form", "<html><head/><body><p>Visit <a href=\"http://femtolab.au.dk\"><span style=\" text-decoration: underline; color:#0057ae;\">femtolab.au.dk</span></a></p></body></html>"))
        self.Ok.setText(_translate("Form", "Close"))

