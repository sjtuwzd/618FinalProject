from alpha_vantage.timeseries import TimeSeries
# from alpha_vantage.techindicators import TechIndicators
# from matplotlib.pyplot import figure
# import matplotlib.pyplot as plt
import json

# Your key here
key = 'ES1O2L3E53BH8P7K'
# Chose your output format, or default to JSON (python dict)
ts = TimeSeries(key)
# ti = TechIndicators(key)

# Get the data, returns a tuple
# aapl_data is a pandas dataframe, aapl_meta_data is a dict
aapl_data, aapl_meta_data = ts.get_daily(symbol='AAPL')
# aapl_sma is a dict, aapl_meta_sma also a dict
# aapl_sma, aapl_meta_sma = ti.get_sma(symbol='AAPL')


# Visualization
# figure(num=None, figsize=(15, 6), dpi=80, facecolor='w', edgecolor='k')
# aapl_data['4. close'].plot()
# plt.tight_layout()
# plt.grid()
# plt.show()

# save
with open('data.json', 'w') as f:
    json.dump(aapl_data, f)
