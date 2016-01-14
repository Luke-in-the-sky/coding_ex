
# coding: utf-8

# # Suggest an option strategy
# 
# We want to buy/sell an option in the SPY market. 
# A good stratey would be one where
#  - Expiration date is some 40 days from now
#  - Probaility of profit and expected returns are high
# 
# To determine our options:
#  1. Load current option chain for SPY
#  2. expDate = next expiration date 40 days from now
#  3. Compute Probability of Profit (POP)
#      - Simulate SPY price evolution untill expDate (Monte Carlo)
#      - Assess strategy: how often have you made a profit?
#      - Check: if you plot Bid vs POP you should see some positive 
#      correlation (if you buy the option, you'll have to pay more 
#      for an option more likely to give you a profit)
#  4. Select best options based on POP and Bid price
# 
# Assuming:
#  - European-type options
#  - strategy: only naked call/put
#  - we are buying the option (i.e. using 'Bid' and not 'Ask')
#  
# To add:
#  - [] Monte Carlo should look at historic data to measure the historic returns of the underlying
#  - [] historical data: should automatically download fresh copy (if older than a week)
#  - [] option chain: after downloading, if all Bids are zero (i.e. not trading hours) ask to keep older data (if available)
#  - [x] implement possibility to *sell* the option (now it's only considering the buy side)
#  - [x] expected return = median(returns)
#  
#  
# Probelms:
#   - [x] Check it pulls reat time data (when US stock exchange is open)
#   - [] I don't like the equation for the Monte Carlo process. Why should we subtract the sigma term from \mu?



import os.path
import math
import pandas
from pandas.io.data import SymbolWarning, RemoteDataError, Options

import numpy 
import datetime
from scipy.stats import norm

my_local_CSV_file = 'C:\\Users\\Luke-Lenovo\\Documents\\GitHub\\coding_ex\\Python - Options strategy\\SPY_YahooFinance.csv'

def side_by_side(*objs, **kwds):
    from pandas.core.common import adjoin
    space = kwds.get('space', 4)
    reprs = [repr(obj).split('\n') for obj in objs]
    print (adjoin(space, *reprs))
    

trading_days_in_a_year = 252


# # Load historical data
# 
# Compute historical volatility and average returns



def load_from_csv(filename = 'SPYopenDayly.csv', verbose=0):
    """ 
    Loading the data file from .csv into dataframe 
    Returns a dataframe object with the read data
    """
    
    # import data in dataframe object
    data = pandas.read_csv(filename)
    data = data.set_index(['Date'])
    if verbose:
        print("Checking keys:")
        print(data.keys())
    
    return data




def compute_hist_vol(full_data, days_in_history_to_consider=40):
    """ 
    Computing historical volatility from hystorical values 
    (CLOSE-to-CLOSE, no dividend, annualized)
    """
    # check you have at least 3 days to work with
    if days_in_history_to_consider < 3:
        print(' --> compute_hist_vol(): need at least 3 days to work with')
        return None

    # define the math for the historical volatility
    try:
        today_closing = full_data.ix[:days_in_history_to_consider]['Close']
        previous_closing = today_closing.shift(-1)
        day_return = (today_closing)/previous_closing
        log_return = day_return.map(lambda x: math.log(x))
    except KeyError:
        print(' --> compute_hist_vol(): could not find the closing prices         in the data')
        return None
 
    # return hist vol
    hist_vol = (log_return.std()/(math.sqrt(1/trading_days_in_a_year)))
    return [hist_vol, day_return.map(lambda x: x-1)]




# ## Retrieve option chain from Yahoo Finance

def retrieve_option_chain(ticker = 'SPY'):
    filename = str(datetime.date.today()) + '_' + ticker + 'optionChain.csv'
    spy_data = pandas.DataFrame()
    if not (os.path.isfile(filename)):
        print (" --> Today's option chain not on system: \
        downloading it now... ")
        try:
            spy_options = Options(ticker, 'yahoo')
            spy_data = spy_options.get_all_data()
            spy_data = spy_data.reset_index()
            spy_data = spy_data.set_index('Symbol')
            spy_data.to_csv(filename, date_format='%Y-%m-%d')
        except:
            print ('---> Error with retrieving data. Please double-check the ticker symbol')
            return spy_data

    spy_data = pandas.read_csv(filename)
    spy_data = spy_data.set_index('Symbol')
    
    return spy_data 

# TODO: delete old Option Chains
#os.remove()




# ## Find an interesting expiry date



import datetime

def next_available_expiration_date(underlyingDataFrame,
                                   how_many_days_from_now = 40):
    '''
    Looks up the given pandas.io.data.Options object and returns 
    the next expiration date (string)
    '''
    # get all unique values in 'Expiry'
    expDates_str = pandas.Series(underlyingDataFrame.Expiry.ravel()).unique()
    # convert values to date obj
    expDates = [ datetime.datetime.strptime(exp_str, '%Y-%m-%d').date()
                for exp_str in expDates_str ] 
    
    expDates.sort()
    for date in expDates:
        if date > date.today() + datetime.timedelta(how_many_days_from_now):
            expiry = date
            break
    return (expiry)




# ## Compute probability of Profit

# #### MonteCarlo simulation
# https://en.wikipedia.org/wiki/Monte_Carlo_methods_in_finance#Monte_Carlo_methods
# 
# To sample a path following the Black-Scholes model from time 0 to T, we chop the time interval into M units of length $\delta t$, and approximate the Brownian motion over the interval $\delta t$ by a single normal variable of mean 0 and variance $\delta t$. This leads to a sample path of
# 
# $$
# S( k\delta t) = S(0) \exp\left( \sum_{i=1}^{k} \left[\left(\mu - \frac{\sigma^2}{2}\right)\delta t + \sigma\varepsilon_i\sqrt{\delta t}\right] \right)
# $$
# 
# for each k between 1 and M. Here each $\varepsilon_i$ is a draw from a standard normal distribution



def evolve_underlying_price(underlying_price = 200, 
                            averageUnderlyingReturns = .0907,
                            volatility = .2, days_to_expiration = 7, 
                            nTrials = 1000):
    '''
    Monte Carlo simulation of stock price evolution.
    Assuming a lognormal distribution of prices, computes 
    price paths and returns a DataFrame(days_to_expir X nTrials):
    each column is the day-by-day evolution of the stock price
    '''
    t = 1/252
    moves = []
    mu = (averageUnderlyingReturns - pow(volatility,2)/2)
    for trial in range(nTrials):
        randArray = numpy.random.normal(0, 1, days_to_expiration-1)
        prices = [underlying_price]
        for eps in randArray:
            change = math.exp(mu*t + volatility*eps*pow(t,0.5))
            prices.append(numpy.round(prices[-1]*change, decimals=2))
        moves.append(prices)

    df = pandas.DataFrame(moves)
    df = df.transpose()
    return df




# #### Compute profit/loss at expiration



def profit_at_expiration_for_single_option(underlyingPriceAtExpiration=100, 
                                           strike=90, optionType='call', 
                                           Price=2):
    underlyingPrice = underlyingPriceAtExpiration
    '''
    Compute the returns from a naked call/put option at expiration given: 
     - the debit/credit payed/gained for it (i.e. 'Price')
     - the underlying price at expiration
     
    Use parameter 'Price' to differenciate between having bought (Price>0) 
    or sold (Price<0)the option.
    '''
    if Price == 0:
        print(' --> profit_at_expiration_for_single_option():         can not work with a Bid price == 0 $')
    optionBought = (Price>0)
    Price = math.fabs(Price)
    if optionType == 'call':
        if optionBought: 
            output = max(underlyingPrice - strike, 0) - Price
        else:
            output = Price -1* max(underlyingPrice - strike, 0)
    elif optionType == 'put':
        if optionBought: 
            output = max(strike - underlyingPrice, 0) - Price
        else:
            output = Price -1* max(strike - underlyingPrice, 0)
    else:
        print(' --> profit_at_expiration_for_single_option():         do not recognize option type')
        output = None
    return output






def strategy_assessment(strike, optionType, price, underlPathsDataFrame, 
                        getOutcomes=False):
    '''    
    Compute the returns from a naked call/put option at expiration given: 
     - the debit/credit payed/gained for it (i.e. 'Price')
     - a simulated paths evolution of its underlying. 
    for each path in the simulation
    
    Returns [POP, list_of_returns]
    '''
    outcome = []
    count = 0
    nTrials = underlPathsDataFrame.shape[1] # i.e. num of simulated paths
    for trial in range(nTrials):
        #(underlyingPriceAtExpiration=100, strike=90, optionType='call', Price=2):
        out = profit_at_expiration_for_single_option(underlPathsDataFrame[trial].iloc[-1],
                                                     strike, optionType, price)
        if getOutcomes:
            outcome.append(out)
        if out > 0:
            count += 1
    
    if getOutcomes:
        #print ('Trials:', nTrials, '   outcomes:', len(outcome))
        return [count / nTrials, outcome]
    else:
        return (count / nTrials)




# ### Compute POP for all (interesting) entries in the option chain



def strategy_buy_naked_callput(underlOptionChain, expirationDate, underlSimulPaths):
    expiry = expirationDate
    selection = underlOptionChain[(underlOptionChain['Expiry']==str(expiry)) & 
                                      (underlOptionChain['Bid']>0)][['Strike', 'Type',
                                                                     'Bid', 'Ask', 'IV']]
    pops = pandas.DataFrame()
    if selection.empty:
        print (' --> Dataframe empty!')
        print ('     Nothing found for', expiry, 'and Ask>0')
    else:
        selection = selection.reset_index()

        for i in range(len(selection)):
            sym, s, t, b, a, iv = selection.loc[i]
            price = b
            pop, returns = strategy_assessment(s,t,price, underlSimulPaths, True)
            median_returns = numpy.median(numpy.asarray(returns))
            pops = pops.append(pandas.Series([sym, s, t, price, pop, median_returns]), ignore_index=True)
        pops.columns = ['Symbol','Strike','Type', 'Price', 'pop', 'mReturns']
    return pops
    




def strategy_sell_naked_callput(underlOptionChain, expirationDate, underlSimulPaths):
    expiry = expirationDate
    selection = underlOptionChain[(underlOptionChain['Expiry']==str(expiry)) & 
                                      (underlOptionChain['Ask']>0)][['Strike', 'Type',
                                                                     'Bid', 'Ask', 'IV']]
    pops = pandas.DataFrame()
    if selection.empty:
        print (' --> Dataframe empty!')
        print ('     Nothing found for', expiry, 'and Ask>0')
    else:
        selection = selection.reset_index()
        for i in range(len(selection)):
            sym, s, t, b, a, iv = selection.loc[i]
            price = -a
            pop, returns = strategy_assessment(s,t,price, underlSimulPaths, True)
            median_returns = numpy.median(numpy.asarray(returns))
            pops = pops.append(pandas.Series([sym, s, t, price, pop, median_returns]), ignore_index=True)
        pops.columns = ['Symbol','Strike','Type', 'Price', 'pop', 'mReturns']
    return pops





