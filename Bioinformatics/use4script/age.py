# -*- coding: utf-8 -*-

import pandas as pd

df = pd.read_excel('Age.xlsx')


def f(x):
    days = 0
    x1 = x.split('岁')
    try:
        if len(x1) > 1:
            year = int(x1[0])
            x = x1[1]
        else:
            year = 0
            x = x1[0]

        x1 = x.split('月')
        if len(x1) > 1:
            month = int(x1[0])
            x = x1[1]
        else:
            month = 0
            x = x1[0]

        x1 = x.split('天')
        if len(x1) > 1:
            day = int(x1[0])
        else:
            day = 0

        days = 365*year + 30*month + day

        return days
    except:
        print(x)
        return 'error'
    
df['days'] = df['年龄'].map(f)

df.to_excel('age.day.xlsx', index=False)