import cPickle as cpkl
import sys

if len(sys.argv) < 2:
    print "Usage: %s pickle file name.\n" % (sys.argv[0])
else:
    with open('data.pkl', 'rb') as f:
        data = cpkl.load(f)

    print data
