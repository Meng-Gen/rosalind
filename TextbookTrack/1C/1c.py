import re
import sys

pattern = sys.stdin.readline().strip()
text = sys.stdin.readline().strip()
print(' '.join([str(m.start()) for m in re.finditer('(?=%s)' % pattern, text)]))
