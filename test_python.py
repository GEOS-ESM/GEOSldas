import datetime
from dateutil.relativedelta import relativedelta

start = datetime.date(2004,1,1)
end = datetime.date(2007,1,1)

date = start + relativedelta(years=1)
delta = date - start
years = delta.days
print(years)
