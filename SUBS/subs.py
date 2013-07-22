import re
import sys

s = sys.stdin.readline().strip()
t = sys.stdin.readline().strip()
print(' '.join([str(m.start()+1) for m in re.finditer('(?=%s)' % t, s)]))
