import timeit
import os
import subprocess

ip = 'python3 sequenceAlign.py'.encode('utf-8')
result = subprocess.run(['/usr/bin/time', 'python3', 'sequenceAlign.py'], stdout=subprocess.PIPE)
result.stdout.decode('utf-8')
print('------------------------------')
print(result.stdout)
print(result.stdout.decode('utf-8'))