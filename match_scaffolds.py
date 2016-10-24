

import sys

# ground truth
with open(sys.argv[1]) as f:
    a = f.readlines()
a = map(lambda x: x.strip(), a)
a = ",".join(a)
a = a.split(">")[1:]
a = map(lambda y: [x for x in y.split(",") if x != ''][1:], a)
a = map(lambda y: map(lambda z: z.split(":::")[0], y), a)
print a



with open(sys.argv[2]) as f:
    b = f.readlines()
b = map(lambda x: x.strip(), b)
b = ",".join(b)
b = b.split(">")[1:]
b = map(lambda y: [x for x in y.split(",") if x != ''][1:], b)
b = map(lambda y: map(lambda z: z.split(":::")[0], y), b)
b = [x for x in b if len(x) > 1]
print b





