#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function

import sys
import numpy as np
import pandas as pd

import random
import datetime

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt

import PyQt4.QtCore as QtCore
import PyQt4.QtGui as QtGui

import lib_functions as mylib

class Example(QtGui.QMainWindow):
    
    def __init__(self):
        super(Example, self).__init__()
        self.tick_str   = 'SPY'
        self.spy_data   = pd.DataFrame()
        self.SimulPaths = pd.DataFrame()
        self.expiration = None
        self.probabilitiesOfProfit = pd.DataFrame()
        self.initUI()
        
    def initUI(self):
    
        # Menu - File
        exitAction = QtGui.QAction('Exit', self) 
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtGui.qApp.quit)
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(exitAction)
        
        # Menu - Help
        aboutAction = QtGui.QAction('About', self) 
        aboutAction.triggered.connect(self.aboutMessage)
        helpMenu = menubar.addMenu('&Help')
        helpMenu.addAction(aboutAction)
        
        # Widgets
        self.central = QtGui.QWidget(self)
        self.setCentralWidget(self.central)
        
        # Ticker Sym [Line Edit]
        self.lbl1 = QtGui.QLabel('Ticker', self.central)
        self.le = QtGui.QLineEdit(self.central)
        self.le.setText('SPY')
        self.btn = QtGui.QPushButton('Get data', self.central)
        self.btn.clicked.connect(self.getTickerSym)
        
        # Strategy [Combo box]
        self.lbl2 = QtGui.QLabel("Strategy", self.central)
        self.combo = QtGui.QComboBox(self.central)
        self.combo.addItem("Buy naked call/put")
        self.combo.addItem("Sell naked call/put")
        self.combo.activated[str].connect(self.munchTheNumbers) 
        
        # Canvas
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        # Layout
        grid = QtGui.QGridLayout()
        grid.setSpacing(10)
        grid.addWidget(self.lbl1    , 1, 0)
        grid.addWidget(self.le      , 1, 1)
        grid.addWidget(self.btn     , 1, 2)
        grid.addWidget(self.lbl2    , 2, 0)
        grid.addWidget(self.combo   , 2, 1)
        grid.addWidget(self.canvas  , 3, 1)
        #grid.addWidget(self.button  , 3, 2)
        grid.addWidget(self.toolbar, 4, 1)
        
        self.central.setLayout(grid)
        self.setGeometry(300, 200, 600, 500)
        self.setWindowTitle('Input dialog')
        self.show()
        self.statusBar().showMessage('Ready')
    
    def aboutMessage(self):
        QtGui.QMessageBox.information(self, 'Message',
            "Option and underlying's data downloaded from Yahoo!\n\n\
            App developped by Luca: suggestions are most welcomed\n\
            at luca.mail.inbox@gmail.com", 
            buttons = QtGui.QMessageBox.Ok)
    
    def getTickerSym(self):        
        self.tick_str = self.le.text()
        self.get_data()
    
    def get_data(self):
        self.statusBar().showMessage('Loading Option data.. (this might take a minute or two)')
        self.spy_data   = mylib.retrieve_option_chain(self.tick_str)
        if self.spy_data.empty:
            reply = QtGui.QMessageBox.warning(self, 'Message',
            'Error while loading option data. \nPlease double-check the ticker symbol is correct', 
            buttons = QtGui.QMessageBox.Ok)
            return
        self.simulatePaths()
            
    def simulatePaths(self):
        price_now       = self.spy_data['Underlying_Price'].iloc[0]
        self.expiration = mylib.next_available_expiration_date(self.spy_data, 40)
        days_to_exp     = self.expiration - datetime.date.today()

        histData        = mylib.load_from_csv(mylib.my_local_CSV_file)
        hist_vol, r     = mylib.compute_hist_vol(histData, 300)
        avAnnualReturns = r.mean() # <== NOT COUNTING FOR DIVIDENDS
        
        self.statusBar().showMessage('Simulating underlying evolution (Monte Carlo)..')
        self.SimulPaths  = mylib.evolve_underlying_price(underlying_price = price_now, 
                        days_to_expiration=days_to_exp.days, nTrials=500, 
                        averageUnderlyingReturns = avAnnualReturns, volatility = hist_vol)
        
        self.munchTheNumbers()        
        self.statusBar().showMessage('Done.')

    def munchTheNumbers(self):
        # Check you have the data you need
        if (self.spy_data.empty or self.SimulPaths.empty or self.expiration==None):
            self.error_no_data_loaded()
            return
            
        # Select strategy to assess
        self.statusBar().showMessage('Assessing strategy..')
        strategy_name = self.combo.currentText()
        if strategy_name == "Buy naked call/put":
            pops = mylib.strategy_buy_naked_callput(underlOptionChain=self.spy_data, 
                        expirationDate=self.expiration, underlSimulPaths=self.SimulPaths)
        elif strategy_name == "Sell naked call/put":
            pops = mylib.strategy_sell_naked_callput(underlOptionChain=self.spy_data, 
                        expirationDate=self.expiration, underlSimulPaths=self.SimulPaths)
        self.probabilitiesOfProfit = pops
        
        self.plot()
            
    def error_no_data_loaded(self):
        self.statusBar().showMessage('Error: no data')
        reply = QtGui.QMessageBox.question(self, 'Message',
                'No data loaded (yet!). \nShould I load the data for %s?' % self.tick_str, 
                QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.Yes)

        if reply == QtGui.QMessageBox.Yes:
            self.get_data()
            
        
    def plot(self):
        self.statusBar().showMessage('Plotting results..')
        data = [random.random() for i in range(10)]
        self.ax = self.figure.add_subplot(111) # create an axis
        #self.ax.hold(False)          # discards the old graph
        plt.cla() #clears an axis, i.e. the currently active axis in the current figure. It leaves the other axes untouched.
        cPops = self.probabilitiesOfProfit[self.probabilitiesOfProfit['Type'] == 'call']
        cPops.plot(ax = self.ax, kind = 'scatter', y='mReturns', x='pop', color='red', label='call')
        cPops = self.probabilitiesOfProfit[self.probabilitiesOfProfit['Type'] == 'put']
        cPops.plot(ax = self.ax, kind = 'scatter', y='mReturns', x='pop', color='blue', label='put')
        plt.legend()
        plt.ylabel('Median returns')
        plt.xlabel('Probability of Profit')
        self.figure.tight_layout()
        self.canvas.draw()      # refresh canvas
        self.statusBar().showMessage('Done.')
        
def main():    
    app = QtGui.QApplication(sys.argv)
    ex = Example()
    #ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()