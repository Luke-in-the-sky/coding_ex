import sys
import string
import logging

from util import mapper_logfile
logging.basicConfig(filename=mapper_logfile, format='%(message)s',
                    level=logging.INFO, filemode='w')

def mapper():
    '''
    For this exercise, compute the average value of the ENTRIESn_hourly column
    for different weather types. Weather type will be defined based on the
    combination of the columns fog and rain (which are boolean values).
    For example, one output of our reducer would be the average hourly entries
    across all hours when it was raining but not foggy.

    Each line of input will be a row from our final Subway-MTA dataset in csv format.
    You can check out the input csv file and its structure below:
    https://www.dropbox.com/s/meyki2wl9xfa7yk/turnstile_data_master_with_weather.csv
   
    Note that this is a comma-separated file.

    This mapper should PRINT (not return) the weather type as the key (use the
    given helper function to format the weather type correctly) and the number in
    the ENTRIESn_hourly column as the value. They should be separated by a tab.
    For example: 'fog-norain\t12345'
   
    Since you are printing the output of your program, printing a debug
    statement will interfere with the operation of the grader. Instead,
    use the logging module, which we've configured to log to a file printed
    when you click "Test Run". For example:
    logging.info("My debugging message")
    Note that, unlike print, logging.info will take only a single argument.
    So logging.info("my message") will work, but logging.info("my","message") will not.
    '''

    # Takes in variables indicating whether it is foggy and/or rainy and
    # returns a formatted key that you should output.  The variables passed in
    # can be booleans, ints (0 for false and 1 for true) or floats (0.0 for
    # false and 1.0 for true), but the strings '0.0' and '1.0' will not work,
    # so make sure you convert these values to an appropriate type before
    # calling the function.
    def format_key(fog, rain):
        return '{}fog-{}rain'.format(
            '' if fog else 'no',
            '' if rain else 'no'
        )

    fog_ColumnNum = 0
    rain_ColumnNum = 0
    entries_ColumnNum = 0
    Length_of_table = 0
   
    for line in sys.stdin:
        data = line.strip().split(",")
      
        if (fog_ColumnNum == 0):
            fog_ColumnNum = data.index('fog')
            rain_ColumnNum = data.index('rain')
            entries_ColumnNum = data.index('ENTRIESn_hourly')
            Length_of_table = len(data)
          
        elif (len(data) == Length_of_table):
            fog = float(data[fog_ColumnNum])
            rain = float(data[rain_ColumnNum])
            weather = format_key(fog, rain)
            print "%s\t%.1f" % (weather, float(data[entries_ColumnNum]))
            #logging.info("%s\t%.1f" % (weather, float(data[entries_ColumnNum])))


mapper()



def reducer():
    '''
    Given the output of the mapper for this assignment, the reducer should
    print one row per weather type, along with the average value of
    ENTRIESn_hourly for that weather type, separated by a tab. You can assume
    that the input to the reducer will be sorted by weather type, such that all
    entries corresponding to a given weather type will be grouped together.

    In order to compute the average value of ENTRIESn_hourly, you'll need to
    keep track of both the total riders per weather type and the number of
    hours with that weather type. That's why we've initialized the variable
    riders and num_hours below. Feel free to use a different data structure in
    your solution, though.

    An example output row might look like this:
    'fog-norain\t1105.32467557'

    Since you are printing the output of your program, printing a debug
    statement will interfere with the operation of the grader. Instead,
    use the logging module, which we've configured to log to a file printed
    when you click "Test Run". For example:
    logging.info("My debugging message")
    Note that, unlike print, logging.info will take only a single argument.
    So logging.info("my message") will work, but logging.info("my","message") will not.
    '''

    riders = 0      # The number of total riders for this key
    num_hours = 0   # The number of hours with this key
    old_key = None
   
    for line in sys.stdin:
        data = line.strip().split("\t")
      
        if len(data)!=2:
            continue
      
        this_key, count = data
        if old_key and (this_key != old_key):
            riders = riders/(num_hours)
            print "%s\t%.8f" % (old_key, riders)
            riders=0
            num_hours=0
      
        old_key = this_key
        num_hours +=1
        riders += float(count)
      
    if old_key != None:
        riders = riders/(num_hours)
        print "%s\t%.8f" % (old_key, riders)

reducer()