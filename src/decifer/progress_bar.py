import sys, os
import datetime

class ProgressBar:

    def __init__(self, total, length, lock, counter, verbose=False, decimals=1, fill=chr(9608), prefix = 'Progress:', suffix = 'Complete'):
        self.total = total
        self.length = length
        self.decimals = decimals
        self.fill = fill
        self.prefix = prefix
        self.suffix = suffix
        self.lock = lock
        self.counter = counter
        assert lock is not None or counter == 0
        self.verbose = verbose

    def progress(self, advance=True, msg=""):
        flush = sys.stderr.flush
        write = sys.stderr.write
        if advance:
            with self.counter.get_lock():
                self.counter.value += 1
        percent = ("{0:." + str(self.decimals) + "f}").format(100 * (self.counter.value / float(self.total)))
        filledLength = int(self.length * self.counter.value // self.total)
        bar = self.fill * filledLength + '-' * (self.length - filledLength)
        rewind = '\x1b[2K\r'
        result = '%s |%s| %s%% %s' % (self.prefix, bar, percent, self.suffix)
        msg = '[{:%Y-%b-%d %H:%M:%S}]'.format(datetime.datetime.now()) + msg
        if not self.verbose:
            toprint = rewind + result + " [%s]" % (msg)
        else:
            toprint = rewind + msg + "\n" + result
        with self.lock:
            #print("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#")
            #print(toprint.encode('utf-8'))
            #print(write(toprint))
            #print("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#")
            #write(toprint.encode('utf-8'))
            x = write(toprint)
            flush()
            if self.counter.value == self.total:
                write("\n")
                flush()
        return True