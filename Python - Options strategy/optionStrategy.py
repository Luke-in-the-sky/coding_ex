#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function

import sys
import os.path
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
        
        # Initialize default values for key objects
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
        openChainFile = QtGui.QAction('Load option-chain', self)
        openChainFile.triggered.connect(self.load_option_chain_file)
        
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openChainFile)
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
        self.lbl0 = QtGui.QLabel('Ticker', self.central)
        self.le = QtGui.QLineEdit(self.central)
        self.le.setText('SPY')
        self.btn = QtGui.QPushButton('Get data', self.central)
        self.btn.clicked.connect(self.getTickerSym)
        
        # Selected Option Chain
        self.lbl_0 = QtGui.QLabel('Selected Chain', self.central)
        self.lbl_1 = QtGui.QLabel(mylib.option_chain_csv_file, self.central)
        
        # Strategy [Combo box]
        self.lbl2 = QtGui.QLabel("Strategy", self.central)
        self.combo = QtGui.QComboBox(self.central)
        self.combo.addItem("Buy naked call")
        self.combo.addItem("Buy naked put")
        self.combo.addItem("Sell naked call")
        self.combo.addItem("Sell naked put")
        self.combo.activated[str].connect(self.check_ready_to_munch_numbers) 
        
        # Canvas
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        # Selected Option Chain
        self.lb3_0 = QtGui.QLabel('Selected Option', self.central)
        self.lb3_1 = QtGui.QTextEdit('This field will display info on your chosen option:')
        self.lb3_1.append('  1. Enter ticker symbol and hit the "Get Data" button')
        self.lb3_1.append('  2. Select a strategy')
        self.lb3_1.append('  3. Click on a point in the graph to find option\' details')
        
        # Layout
        grid = QtGui.QGridLayout()
        grid.setSpacing(10)
        grid.addWidget(self.lbl0    , 1, 0)
        grid.addWidget(self.le      , 1, 1)
        grid.addWidget(self.btn     , 1, 2)
        grid.addWidget(self.lbl_0   , 2, 0)
        grid.addWidget(self.lbl_1   , 2, 1)        
        grid.addWidget(self.lbl2    , 3, 0)
        grid.addWidget(self.combo   , 3, 1)
        grid.addWidget(self.canvas  , 4, 1)
        grid.addWidget(self.toolbar , 5, 1)
        grid.addWidget(self.lb3_0   , 6, 0)
        grid.addWidget(self.lb3_1   , 6, 1) 
        
        self.central.setLayout(grid)
        self.setGeometry(300, 50, 800, 500)
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
            QtGui.QMessageBox.warning(self, 'Message',
            'Error while donwloading option data. \nPlease double-check the ticker symbol is correct', 
            buttons = QtGui.QMessageBox.Ok)
            return
        self.check_ready_to_munch_numbers()
        
    def ask_to_reset_sim_paths(self):
        reply = QtGui.QMessageBox.question(self, 'Message',
            'Should I reset the simulated paths?', 
            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.Yes)
        if reply == QtGui.QMessageBox.Yes:
            self.SimulPaths = pd.DataFrame()  # Reset the simulated paths
            if hasattr(self, 'ax_sim'):
                self.ax_sim.cla()
                self.canvas.draw()      # refresh canvas
        
    def load_option_chain_file(self):
        self.statusBar().showMessage('Loading Option Chain file')
        self.ask_to_reset_sim_paths()
        csv_filename = QtGui.QFileDialog.getOpenFileName(self, 'Open csv file')
        self.spy_data = mylib.retrieve_option_chain_from_file(csv_filename)
        self.check_ready_to_munch_numbers()
        
    def simulatePaths(self):
        price_now       = self.spy_data['Underlying_Price'].iloc[0]
        self.expiration = mylib.next_available_expiration_date(self.spy_data, 40)
        days_to_exp     = self.expiration - datetime.date.today()
        history_start_date = datetime.date.today() - datetime.timedelta(600)
        try:
            histData = mylib.download_historical_data(self.tick_str, history_start_date)
        except:
            self.error_hist_data_not_found()
            return
            
        hist_vol, r     = mylib.compute_hist_vol(histData, mylib.days_in_history_to_consider)
        avAnnualReturns = r.mean() # <== NOT COUNTING FOR DIVIDENDS
        
        self.statusBar().showMessage('Simulating underlying evolution (Monte Carlo)..')
        self.SimulPaths  = mylib.evolve_underlying_price(underlying_price = price_now, 
                        days_to_expiration=days_to_exp.days, nTrials=500, 
                        averageUnderlyingReturns = avAnnualReturns, volatility = hist_vol) 
        # Plot it
        self.plot_simulated_paths()
    
    
    def check_ready_to_munch_numbers(self):
        self.lbl_1.setText(mylib.option_chain_csv_file)
        all_fine = True
        if self.spy_data.empty:
            all_fine = False
            self.error_no_data_loaded()
            return
        sumBid = (self.spy_data['Bid'] **2).sum()
        sumAsk = (self.spy_data['Ask'] **2).sum()
        if sumBid == 0 and sumAsk == 0:
            all_fine = False
            QtGui.QMessageBox.information(self, 'Message',
                'All Bid and Ask prices are set to zero in this option chain.\
                \nPerhaps the data was downloaded outside of the US stock exchange trading hours.\
                \nTry loading option chains from the File menu', 
                QtGui.QMessageBox.Ok)     
            return
        if self.SimulPaths.empty:
            self.simulatePaths()
        
        self.munchTheNumbers()
  

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
        elif strategy_name == 'Buy naked call':
            pops = mylib.strategy_naked_option(underlOptionChain = self.spy_data, 
                   expirationDate = self.expiration, underlSimulPaths=self.SimulPaths,
                   bidOrAsk='Bid', type='call')
        elif strategy_name == 'Buy naked put':
            pops = mylib.strategy_naked_option(underlOptionChain = self.spy_data, 
                   expirationDate = self.expiration, underlSimulPaths=self.SimulPaths,
                   bidOrAsk='Bid', type='put')
        elif strategy_name == 'Sell naked call':
            pops = mylib.strategy_naked_option(underlOptionChain = self.spy_data, 
                   expirationDate = self.expiration, underlSimulPaths=self.SimulPaths,
                   bidOrAsk='Ask', type='call')
        elif strategy_name == 'Sell naked put':
            pops = mylib.strategy_naked_option(underlOptionChain = self.spy_data, 
                   expirationDate = self.expiration, underlSimulPaths=self.SimulPaths,
                   bidOrAsk='Ask', type='put')
        self.probabilitiesOfProfit = pops
        self.plot()
            
    def error_no_data_loaded(self):
        self.statusBar().showMessage('Error: no option data')
        reply = QtGui.QMessageBox.information(self, 'Message',
                'No data loaded (yet!).', 
                QtGui.QMessageBox.Ok)

    
    def error_hist_data_not_found(self):
        self.statusBar().showMessage('Error: no historical data')
        QtGui.QMessageBox.warning(self, 'No historical data',
            'No historical data found. \nPlease check you have it at:\n %s' % mylib.my_local_CSV_file, 
            QtGui.QMessageBox.Ok)        
            
    def plot_simulated_paths(self):
        self.statusBar().showMessage('Plotting Simulated paths..')
        self.ax_sim = self.figure.add_subplot(121) # create an axis
        plt.cla()
        if self.SimulPaths.empty:
            print('No simulated paths to be plotted')
            return
        self.ax_sim.set_title('Prices paths')
        self.ax_sim.set_ylabel('Underl. price')
        self.ax_sim.set_xlabel('Days')
        self.SimulPaths.plot(ax = self.ax_sim)
        self.ax_sim.legend_ = None
        self.figure.tight_layout()
        self.canvas.draw()      # refresh canvas
                
    
    def plot(self):
        self.statusBar().showMessage('Plotting results..')
        self.ax = self.figure.add_subplot(122) # create an axis
        plt.cla() #clears an axis, i.e. the currently active axis in the current figure. It leaves the other axes untouched.
        if self.probabilitiesOfProfit.empty:
            QtGui.QMessageBox.information(self, 'Message',
                'Nothing to be plotted.', 
                QtGui.QMessageBox.Ok) 
        else:
            Pops = self.probabilitiesOfProfit
            Pops.plot(ax = self.ax, kind = 'scatter', y='mReturns', x='pop', color='blue', picker=True)
            
            def pick_option_from_plot(event):
                ind = event.ind
                data = self.probabilitiesOfProfit
                clicked_pop = (np.take(Pops['pop'], ind)).iloc[0]
                clicked_ret = (np.take(Pops['mReturns'], ind)).iloc[0]
                selected_option = Pops[(Pops['pop'] == clicked_pop) & (Pops['mReturns'] == clicked_ret) ]
                print(selected_option)
                sym, strike, type, price, pop, mRet = selected_option.iloc[0]
                selected_option_str = 'SYM: %s \nPrice:    %.2f $ \nPoP:     %.2f %% \nMedian Returns: %.2f $' % (
                    sym, price, pop*100, mRet)
                self.lb3_1.setText(selected_option_str)
            
            self.ax.set_title('Strategy Assesment')
            self.ax.set_ylabel('Median returns')
            self.ax.set_xlabel('Probability of Profit')
            self.figure.tight_layout()
            self.canvas.draw()      # refresh canvas
            self.figure.canvas.mpl_connect('pick_event', pick_option_from_plot)
            self.statusBar().showMessage('Done.')


def main():    
    app = QtGui.QApplication(sys.argv)
    ex = Example()
    #ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()