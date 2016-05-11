#!/usr/bin/python3

#   Copyright 2016 Anders Aspegren Søndergaard / Femtolab, Aarhus University
#
#   This file is part of Alignment calculator.
#
#   Alignment calculator is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   Alignment calculator is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with Alignment calculator. If not, see <http://www.gnu.org/licenses/>.


import os,sys
import PyQt5
import PyQt5.QtCore
import PyQt5.QtGui
import PyQt5.QtWidgets
from PyQt5.QtWidgets import QWidget, QMessageBox, QApplication, QFileDialog
from PyQt5.QtCore import QThread, QObject
from PyQt5.QtCore import pyqtSlot, pyqtSignal
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import numpy
from _gui import Ui_MainWindow
from _about import Ui_Form as aboutForm
from _precalculate import Ui_precalculate as precalculateForm
import time
import configparser,config,boltzmann
import tempfile
from tempfile import NamedTemporaryFile
import utils
import gc
import multiprocessing
import U2dcalc


#cmd = sys.executable;


class GUI(PyQt5.QtWidgets.QMainWindow):
    def __init__(self,hardcore,scratchdir):
        super().__init__();
        self.hardcore = hardcore;
        self.scratchdir = scratchdir
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self);
        
        self.moleculeConfigFile = "conf/molecules.conf";

        _translate = PyQt5.QtCore.QCoreApplication.translate
        self.ui.moleculesBox.clear();
        self.ui.moleculesBox.addItems(self.list_molecules());
        self.ui.moleculesBox.model().sort(0);
        self.ui.moleculesBox.model().insertRow(0);
        self.ui.moleculesBox.setItemText(0,"Load...");
        self.ui.moleculesBox.setCurrentIndex(0);


        self.ui.mybutn.clicked.connect(self.Go);
        self.ui.saveTrace.triggered.connect(self.saveTrace);
        self.ui.actionPrecalculate.triggered.connect(self.precalculate);
        self.ui.actionClose_figures.triggered.connect(self.close_figs);
        self.ui.actionExit.triggered.connect(self.quit_program);
        self.ui.actionHelp.triggered.connect(self.help);
        self.ui.actionAbout.triggered.connect(self.about);
        self.ui.moleculesBox.activated.connect(self.moleculeSelected);
        self.ui.alpha_par.textChanged.connect(self.calculateDalpha);
        self.ui.alpha_perp.textChanged.connect(self.calculateDalpha);
        
        self.ui.Aconst.textChanged.connect(self.update_num_states);
        self.ui.Bconst.textChanged.connect(self.update_num_states);
        self.ui.Bconst.textChanged.connect(self.update_timestep);
        self.ui.Temperature.textChanged.connect(self.update_num_states);
        self.ui.abundanceEven.textChanged.connect(self.update_num_states);
        self.ui.abundanceOdd.textChanged.connect(self.update_num_states);
        self.ui.percentile.textChanged.connect(self.update_num_states);
        self.ui.forceDT.clicked.connect(self.force_timestep_click);
        self.ui.Jmax.textChanged.connect(self.update_timestep);
        self.ui.cos2d.clicked.connect(self.update_timestep);

        self.ensemble = [];

        self.show();
    
    def update_num_states(self,checked=0):
        A = default("0",self.ui.Aconst.text());
        B = self.ui.Bconst.text();
        T = default("0",self.ui.Temperature.text());
        even = default("1",self.ui.abundanceEven.text());
        odd = default("1",self.ui.abundanceOdd.text());
        Jmax = default("140",self.ui.Jmax.text());
        percentile = default("99.9",self.ui.percentile.text());
        try:
            A = 2*numpy.pi*1e9*float(A);
            B = 2*numpy.pi*1e9*float(B);
            T = float(T);
            even = float(even);
            odd = float(odd);
            Jmax = int(Jmax);
            percentile = float(percentile)/100;
            m = config.molecule("MOL");
            m.A = A;
            m.B = B;
            m.even = even;
            m.odd = odd;
            if (T >= 0 and even >= 0 and odd >= 0 and even+odd != 0 and percentile > 0 and percentile < 1):
                ensemble = boltzmann.thermal_JKM_ensemble(T,m,Jmax,percentile);
                self.ensemble = ensemble;
                num = len(ensemble);
                self.ui.num_ensemble_states.setText(str(num));
            else:
                self.ui.num_ensemble_states.setText("N/A");
                self.ensemble = [];
        except Exception as e:
            self.ui.num_ensemble_states.setText("N/A");
            self.ensemble = [];


    def update_timestep(self,checked=False):
        if (not self.ui.timestep.isEnabled()):
            self.ui.timestep.setText("");
            try:
                cos2d = self.ui.cos2d.isChecked();
                Jmax = max(2,int(default("140",self.ui.Jmax.text())));
                B = float(self.ui.Bconst.text());
                dt = 1000*utils.nice_time_step(B,Jmax,for_cos2d=cos2d);
                if dt >= 0.01 and dt < 1000:
                    self.ui.timestep.setText('%.3f' % dt);
                else:
                    self.ui.timestep.setText('%.2e' % dt);
            except:
                pass;


    def force_timestep_click(self,checked):
        self.ui.timestep.setEnabled(checked)


    def closeEvent(self, event):
        self.close_figs();
        self.calc = None;
        self.last_result = None;
        event.accept();

    def quit_program(self):
        self.close_figs();
        self.calc = None;
        self.last_result = None;
        self.close();

    def close_figs(self):
        plt.close('all');

    def list_molecules(self):
        result = [];
        try:
            if (not os.path.exists(self.moleculeConfigFile)):
                raise RuntimeError("Could not load molecules from the file: "+self.moleculeConfigFile);
            parser = configparser.ConfigParser();
            parser.read(self.moleculeConfigFile);
            result = parser.sections();
        except Exception as e:
            self.noMoleculesWarningMsg = str(e);
            timer = PyQt5.QtCore.QTimer(self);
            timer.timeout.connect(self.noMoleculesWarning);
            timer.setSingleShot(True)
            timer.start(0);

        return result

    def noMoleculesWarning(self):
        QMessageBox.warning(self,'Failed to load predefined molecules',self.noMoleculesWarningMsg);

    def disable_Go_button(self):
        self.ui.mybutn.setEnabled(False);
        self.ui.mybutn.setText('Calculating..');

    def enable_Go_button(self):
        self.ui.mybutn.setEnabled(True);
        self.ui.mybutn.setText('Go!');
    
    def disable_save_option(self):
        self.ui.saveTrace.setEnabled(False);
    def enable_save_option(self):
        self.ui.saveTrace.setEnabled(True);

    def Go(self):
        params = self.validate_input();
        if (params):
            self.disable_Go_button();
            self.disable_save_option();

            self.calc = calculatron(params,self.scratchdir);
            # In order to send signals back and forth in a thread safe manner,
            # we must use QThreads.
            self.calc.finished.connect(self.calcDone);
            self.calc.start();

    def validation_error(self,message):
        if (self.hardcore):
            QMessageBox.critical(self,'Validation error',"Moron!");
        else:
            QMessageBox.critical(self,'Validation error',message);
        raise ValueError(message);
    
    def validator_to_float(self,field,num):
        try:
            val = float(num);
            if (val<0):
                raise ValueError("");
            return val;
        except ValueError:
            self.validation_error(field + " must be a nonnegative number.");
            raise
 
    def validator_to_int(self,field,num):
        try:
            val = int(num);
            if (val < 0):
                raise ValueError("");
            return val;
        except ValueError:
            self.validation_error(field + " must be a nonnegative integer.");
            raise
            
    def validator_to_floats(self,field,nums,n_floats=0):
        try:
            vals = [float(num) for num in nums.split(",")]
            if (any([val<0 for val in vals])):
                raise ValueError("");
            if (len(vals) == 1 and n_floats > 1):
                vals *= n_floats; # Make sure we have n_floats floats
            if (len(vals) != n_floats and n_floats > 0):
                raise RuntimeError("");
            return vals;
        except ValueError:
            self.validation_error(field + " must be a nonnegative number.");
            raise
        except RuntimeError:
            if (n_floats > 1):
                self.validation_error(field + " must contain one or "+ str(n_floats) + " comma separated nonnegative numbers.");
            else:
                self.validation_error(field + " must contain one nonnegative number.");
            raise ValueError()

    def validate_input(self):
        try: 
            p = dict()
            Aconst = default("0",self.ui.Aconst.text());
            Bconst = self.ui.Bconst.text()
            alpha_par = self.ui.alpha_par.text()
            alpha_perp = self.ui.alpha_perp.text()
            FWHM = default("300",self.ui.pulseDuration.text())
            I0 = default("10",self.ui.pulseIntensity.text())
            pumpWaist = default("35",self.ui.pumpWaist.text())
            probeWaist = default("25",self.ui.probeWaist.text())
            Nshells = default("1",self.ui.Nshells.text());
            if (Nshells == "0"):
                self.ui.Nshells.setText("1");
                Nshells = "1";
            t0 = default("0",self.ui.t0.text());
            Jmax = default("140",self.ui.Jmax.text());
            temperature = default("0",self.ui.Temperature.text());
            abundanceEven = default("1",self.ui.abundanceEven.text());
            abundanceOdd = default("1",self.ui.abundanceOdd.text());
            percentile = default("99.9",self.ui.percentile.text());
            anisotropy = default("1.0",self.ui.ELfactor.text());
            J = default("0",self.ui.J.text());
            K = default("0",self.ui.K.text());
            M = default("0",self.ui.M.text());
            cos2d = self.ui.cos2d.isChecked();

            p["Aconst"] = self.validator_to_float("A constant",Aconst)
            p["Bconst"] = self.validator_to_float("B constant", Bconst)
            p["a_par"] = self.validator_to_float("Parallel polarizability", alpha_par)
            p["a_perp"] = self.validator_to_float("Perpendicular polarizability", alpha_perp)

            p["t0"] = [i*1e-12 for i in self.validator_to_floats("T0",t0)]
            num_pulses = len(p["t0"]);
    
            p["I_0"] = [i*1e16 for i in self.validator_to_floats("Peak intensity", I0,num_pulses)];
            p["FWHM"] = [i*1e-15 for i in self.validator_to_floats("Pulse duration", FWHM, num_pulses)];

            p["pumpWaist"] = [i*1e-6 for i in self.validator_to_floats("Pump waist",pumpWaist,num_pulses)];

            # Note: only pump waist must be converted to SI units.
            # The calculation program takes the probe waist in micro meters.
            p["probeWaist"] = self.validator_to_float("Probe waist",probeWaist);
            p["Nshells"] = self.validator_to_int("#FVA shells",Nshells);
            p["Jmax"] = self.validator_to_int("Jmax",Jmax);
    
            if (p["pumpWaist"] == 0 or p["probeWaist"] == 0):
                self.validation_error("Waist size must be larger than 0.");
            
            initial_state = self.ui.InitConditionsTab.currentIndex();
            p["Boltzmann"] = False;
            p["singleState"] = False;
            if (initial_state == 0):
                p["Boltzmann"] = True;
                p["temperature"] = self.validator_to_float("Temperature",temperature);
                p["abundanceEven"] = self.validator_to_float("Even abundance",abundanceEven);
                p["abundanceOdd"] = self.validator_to_float("Even abundance",abundanceOdd);
                p["percentile"] = self.validator_to_float("Percentile", percentile)/100;
                if (p["percentile"] <= 0 or p["percentile"] >= 1):
                    self.validation_error("Percentile must be between 0 and 100 %, exclusive.");
                p["anisotropy"] = self.validator_to_float("E-L anisotropy", anisotropy);
            elif (initial_state == 1):
                p["singleState"] = True;
                p["J"] = self.validator_to_int("J",J);
                try:
                    p["K"] = int(K);
                    p["M"] = int(M);
                except ValueError:
                    self.validation_error("K and M must be integers.");
                if (p["Aconst"] == 0 and p["K"] != 0):
                    self.validation_error("Linear molecules have K=0");
                J = p["J"];
                K = p["K"];
                M = p["M"];
                if (J<0 or abs(M)>J or abs(K)>J):
                    self.validation_error("Triangle inequality violated.");
            else:
                self.validation_error("Invalid initial state tab");
    
            p["cos2d"] = cos2d;
    
            if (self.ui.forceDT.isChecked()):
                dt = default("0",self.ui.timestep.text());
                dt = self.validator_to_float("Time step",dt);
                p["dt"] = dt;



            return p;
        except ValueError:
            return False;


    def calcDone(self):
        self.enable_Go_button();
        if (self.calc.error):
            QMessageBox.critical(self,'Calculation error',str(self.calc.error));
        else:
            self.last_result = self.calc.result;
            self.enable_save_option();
            self.present_results();
    
    
    def present_results(self):
        npzfile = self.last_result;
        t = npzfile["t"];
        cos2 = npzfile["cos2"];
        cos2d = npzfile["cos2d"];
        fig = plt.figure();
        plt.plot(t*1e12,cos2);
        plt.xlabel('Time [ps]');
        plt.ylabel('<cos^2 theta>');
        plt.ion();
        fig.show();
        if (len(cos2d) == len(cos2)):
            fig = plt.figure();
            plt.plot(t*1e12,cos2d);
            plt.xlabel('Time [ps]');
            plt.ylabel('<cos^2 theta 2D>');
            plt.ion();
            fig.show();
        if ('Javg' in npzfile.keys()):
            fig = plt.figure();
            Javg = npzfile["Javg"];
            std = npzfile["std"];
            tp = t[:len(Javg)]*1e12;
            percentile_999 = npzfile["percentile_999"];
            plt.plot(tp,Javg,'k-',label='<J>');
            plt.plot(tp,Javg+std,'k--',label='<J>+std(J)');
            plt.plot(tp,Javg-std,'k--',label='<J>-std(J)');
            plt.plot(tp,percentile_999,'r-.',label='99.9 percentile');
            plt.legend(loc=0);
            plt.xlabel('Time [ps]')
            plt.ylabel('J');
            plt.ion();
            fig.show();


    def saveTrace(self,checked):
        #diag = QFileDialog.getSaveFileName(self, "Select destination", "./", "Comma Separated Values (*.csv)");
        diag = QFileDialog(self);
        diag.setAcceptMode(QFileDialog.AcceptSave) #Save file, not open one
        diag.setNameFilter("Comma Separated Values (*.csv);;Space separated Values (*.csv)");
        diag.setDefaultSuffix("csv"); # Make sure selected files end in .csv
        diag.exec();
        try:
            filename = diag.selectedFiles()[0];
        except IndexError:
            filename = '';
        user_filter = diag.selectedNameFilter();
        if (user_filter == "Space separated Values (*.csv)"):
            delimiter = " ";
        else:
            delimiter = ",";

        if (filename != '' and not os.path.isdir(filename)):

            npzfile = self.last_result;
            t = npzfile["t"];
            cos2 = npzfile["cos2"];
            cos2d = npzfile["cos2d"];
                
            extra_header = [];
            extra_columns = [];
            if ('Javg' in npzfile.keys()):
                Javg = npzfile["Javg"];
                std = npzfile["std"];
                percentile_999 = npzfile["percentile_999"];
                extra_header = ["<J>","std(J)","J_99.9%"];
                extra_columns = [Javg,std,percentile_999];

            utils.save_to_csv(filename,t,cos2,cos2d,extra_header,extra_columns,delimiter);

    
    def moleculeSelected(self,index):
        if (index != 0):
            name = self.ui.moleculesBox.itemText(index);
            molecule = config.molecule(name,self.moleculeConfigFile);
            self.ui.Aconst.setText(str(molecule.A_raw));
            self.ui.Bconst.setText(str(molecule.B_raw));
            self.ui.alpha_par.setText(str(molecule.alpha_par_volume));
            self.ui.alpha_perp.setText(str(molecule.alpha_perp_volume));
            self.ui.abundanceEven.setText(str(molecule.even));
            self.ui.abundanceOdd.setText(str(molecule.odd));
        else:
            self.ui.Aconst.setText("");
            self.ui.Bconst.setText("");
            self.ui.alpha_par.setText("");
            self.ui.alpha_perp.setText("");
            self.ui.abundanceEven.setText("");
            self.ui.abundanceOdd.setText("");
 
            

    def calculateDalpha(self):

        try:
            a_par = self.ui.alpha_par.text();
            a_perp = self.ui.alpha_perp.text();
            if (a_par == ""):
                a_par = "0";
            if (a_perp == ""):
                a_perp = "0";
            
            
            self.ui.deltaAlpha.setText(str(round(float(a_par)-float(a_perp),5)));
        except ValueError:
            self.ui.deltaAlpha.setText("");

#        try:
#            result = QMessageBox.information(QWidget(),'Information','Some informative text');
#            print(result);
#            result = QMessageBox.warning(QWidget(),'Warn','WARNONNN!');
#            result = QMessageBox.warning(QWidget(),'Warn','WARNONNN!');
#            print(result);
#            self.enable_save_option();
#        finally:
#            self.enable_Go_button();
#

    def about(self):
        #self.about_w = aboutWidget(self)
        q = aboutWidget(self)

    def precalculate(self):
        ensemble = self.ensemble;
        Jmax = self.ui.Jmax.text();
        try:
            Jmax = max(0,int(Jmax));
        except:
            Jmax = 140;

        Kmax = max([0] + [m[2] for m in ensemble]);
        Mmax = max([0] + [m[3] for m in ensemble]);
            
        q = precalculateWidget(self,Jmax,Kmax,Mmax)

    def help(self):
        helpmsg = "It may be helpful to hover the mouse cursor over\n";
        helpmsg += "items and text labels that you don't understand.";
        QMessageBox.information(self,'Nonadiabatic help',helpmsg);

def default(the_default_value,the_input):
    if (the_input == ""):
        return the_default_value;
    else:
        return the_input;

class aboutWidget(PyQt5.QtWidgets.QDialog):
    def __init__(self,MainWindow):
        super().__init__(MainWindow);
        self.ui = aboutForm();
        self.ui.setupUi(self)
        self.ui.Ok.clicked.connect(self.close);
        self.show();

class precalculateWidget(PyQt5.QtWidgets.QDialog):
    def __init__(self,MainWindow,Jmax,Kmax,Mmax):
        super().__init__(MainWindow);
        self.ui = precalculateForm();
        self.ui.setupUi(self)

        self.ui.cancelButton.clicked.connect(self.close);
        self.ui.calculateButton.clicked.connect(self.calculate);
        self.ui.deleteButton.clicked.connect(self.delete);

        self.ui.Jmax.setText(str(Jmax));
        self.ui.Kmax.setText(str(Kmax));
        self.ui.Mmax.setText(str(Mmax));

        self.ui.Jmax.textChanged.connect(self.update_size_req);
        self.ui.Kmax.textChanged.connect(self.update_size_req);
        self.ui.Mmax.textChanged.connect(self.update_size_req);
        
        self.update_size_req();
        self.update_available();
        self.show();

    def update_available(self):
        size = U2dcalc.cache_size();
        if (size is None):
            self.ui.deleteButton.setEnabled(False);
            self.ui.availableLabel.setText("None.");
        else:
            Jmax,Kmax,Mmax = size;
            self.ui.deleteButton.setEnabled(True);
            self.ui.availableLabel.setText("Jmax = {}, Kmax = {}, Mmax = {}.".format(Jmax,Kmax,Mmax));


    def calculate(self):
        inputs = self.validate_input();
        if (inputs is not False):
            Jmax,Kmax,Mmax = inputs;
            self.setEnabled(False);
            self.update();
            QMessageBox.information(self,'Responsiveness','This user interface will be irresponsive while the calculation is carried out.\n\nSorry about that!');
            try:
                U2dcalc.set_num_threads(multiprocessing.cpu_count());
                U2dcalc.precalculate_matrix_elements(Jmax,Kmax,Mmax);
                self.update_available();
                QMessageBox.information(self,'Success!','Calculation done.');
                self.close();
            except BaseException as e:
                QMessageBox.critical(self,'Failed to calculate matrix elements',str(e));
                self.setEnabled(True);
        else:
            QMessageBox.critical(self,'Validation error',"Invalid input");

    def delete(self):
        U2dcalc.drop_cache();
        self.update_available();

    def validate_input(self):
        inputs = [];
        inputs.append(self.ui.Jmax.text());
        inputs.append(self.ui.Kmax.text());
        inputs.append(self.ui.Mmax.text());

        try:
            inputs = [int(i) if i != "" else 0 for i in inputs];
            for i in inputs:
                if (i < 0):
                    return False;
            return inputs;
        except:
            return False;
    
    def update_size_req(self):
        inputs = self.validate_input();
        if (inputs):
            Jmax,Kmax,Mmax = inputs;
            size = round(8*2*(Jmax+1)**2*(Kmax+1)*(Mmax+1)/1024**2);
            self.ui.sizeLabel.setText(str(size));
            self.ui.calculateButton.setEnabled(True);
        else:
            self.ui.calculateButton.setEnabled(False);




class calculatron(QThread):
    def __init__(self,params,scratchdir):
        super().__init__();
        self.result = None;
        self.error = None;
        self.params = params;
        self.scratchdir = scratchdir
        self.cmd = sys.executable;
        if (' ' in self.cmd):
            self.cmd = '"' + self.cmd + '"';

        p = self.params;
        for name in ["t0","I_0", "FWHM","pumpWaist"]:
            p[name] = ', '.join(str(i) for i in p[name]);

    def run(self):
        #print(self.params)
        try:
            with NamedTemporaryFile(dir=self.scratchdir,delete=False) as conffile:
                self.create_config_file(conffile);
                conffile_name = conffile.name;

            if (self.params["Boltzmann"]):
                self.result = self.calculateBoltzmann(conffile_name);
            elif (self.params["singleState"]):
                self.result = self.calculateSingleState(conffile_name);
            else:
                self.error = "Don't know what to do, giving up! :-(";
        except Exception as e:
            self.error = str(e);
        finally:
            pass;

    def create_config_file(self,conffile):
        p = self.params;
        if (not self.params["Boltzmann"]):
            p["abundanceOdd"] = 1;  # These are not actually used when
            p["abundanceEven"] = 1; # not doing boltzmann averaging
        conffile_text  = "[MOL]\n"
        conffile_text += "A = {}\nB = {}\n".format(p["Aconst"],p["Bconst"])
        conffile_text += "alpha_par_volume = {}\n".format(p["a_par"])
        conffile_text += "alpha_perp_volume = {}\n".format(p["a_perp"])
        conffile_text += "odd = {}\neven = {}\n\n".format(p["abundanceOdd"],p["abundanceEven"])

        conffile_text += "[pulses]\n"
        conffile_text += "I_max = {}\nFWHM = {}\n".format(p["I_0"],p["FWHM"])
        conffile_text += "t = {}\nwaist = {}\n".format(p["t0"],p["pumpWaist"])
        
        conffile.file.write(conffile_text.encode('utf8'))
        conffile.file.flush()

    def calculateBoltzmann(self,conffile_name):
        p = self.params;
        
        cmd = self.cmd;

        args = [cmd, "./cos2_thermal.py MOL", conffile_name, conffile_name]
        args += ["-T", p["temperature"], "--Jmax", p["Jmax"]]
        args += ["--Nshells", p["Nshells"], "--probe_waist", p["probeWaist"]]
        args += ["--percentile", p["percentile"]];
        args += ["--anisotropy", p["anisotropy"]];
        if (p["cos2d"]):
            args += ["--cos2d"]
 
        return self.alignment_trace_command(args);

    def calculateSingleState(self,conffile_name):
        p = self.params;

        cmd = self.cmd;

        args  = [cmd, "./cos2_calc.py MOL", p["J"], p["K"], p["M"]]
        args += [conffile_name,conffile_name, "--Jmax",p["Jmax"]]
        args += ["--Nshells", p["Nshells"], "--probe_waist", p["probeWaist"]]
        if (p["cos2d"]):
            args += ["--cos2d"]

        return self.alignment_trace_command(args);

    def alignment_trace_command(self,args):
        with NamedTemporaryFile(suffix=".npz",dir=self.scratchdir) as resultfile:
            resultfile_name = resultfile.name;

        args += ["--filename",resultfile_name]
        try:
            args += ["--dt", str(self.params["dt"])];
        except:
            pass;
        command = ' '.join(str(i) for i in args);
        print("Executing command:")
        print("");
        print(command);
        print("");
        if (os.system(command) != 0):
            raise RuntimeError("Failed to calculate alignment trace. See console.");
        res = numpy.load(resultfile.name);
        return res;

 
if __name__ == '__main__':
    status = 0;
    with tempfile.TemporaryDirectory() as tmpdirname:
        hardcore = False;
        #hardcore = True
        application = PyQt5.QtWidgets.QApplication(sys.argv);
        gui = GUI(hardcore,tmpdirname);
        status = application.exec_();
        # On windows, for some reason, you are not allowed to delete a file
        # that is open. Orphan all references to open files in the hope that
        # the garbage collector will save our day:
        del gui;
        gc.collect();
    sys.exit(status);

