''' This code obtains the running time from a
gaussian.log file.

python time-gaussian.py arg1

arg1: log file of the completed run.

example:
python time-gaussian.py optimization.log
'''

import subprocess
import sys


def time(keyword):
    p = subprocess.Popen("grep '"+keyword+"' "+sys.argv[1],
                         stdout=subprocess.PIPE,
                         shell=True)
    out = p.communicate()[0].decode('ascii')
    out = out.split()
    start = out.index('at')

    month = out[start+2]
    day = int(out[start+3])
    time = out[start+4].split(':')
    hour = int(time[0])
    minu = int(time[1])
    seco = int(time[2])

    return month, day, hour, minu, seco


t_i = time('Leave Link    1')
t_f = time('Normal termination of Gaussian')

if t_i[0] == t_f[0]:
    days = (t_f[1]-t_i[1])*24*3600
    hours = (t_f[2]-t_i[2])*3600
    minus = (t_f[3]-t_i[3])*60
    secos = t_f[4]-t_i[4]
    total = days + hours + minus + secos

    print("Time in seconds= ", total)
    print("Time in minutes= ", total/60)
    print("Time in hours= ", total/3600)

else:
    print('sorry, I cannot help you, modify me to compute \
           changes of months')
