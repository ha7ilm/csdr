#!/usr/bin/python

import os, time, signal
from subprocess import *
#https://bugs.python.org/issue1652

def p(x):
    global printcmds
    if printcmds: print x
    return check_output(x, shell=True)

printcmds=True


def genfiles(snr):
    cmd="""(while true; do echo -n 'CQ CQ CQ DE HA7ILM HA7ILM HA7ILM PSE K '; done) | \
csdr psk31_varicode_encoder_u8_u8 | \
tee /s/bpsk31_testin | \
csdr differential_encoder_u8_u8 | \
csdr psk_modulator_u8_c 2 | \
csdr psk31_interpolate_sine_cc 256 | \
csdr awgn_cc %d | \
csdr timing_recovery_cc GARDNER 256 0.5 2 --add_q | \
csdr dbpsk_decoder_c_u8 | \
dd bs=1024 count=10 of=/s/bpsk31_testout
"""%snr
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    if printcmds: print cmd
    os.system(cmd)

def getminsize():
    return min(os.path.getsize("/s/bpsk31_testout"), os.path.getsize("/s/bpsk31_testin"))

def mkdiff(shift):
    if shift==0:
        return int(p("cmp -l /s/bpsk31_testin /s/bpsk31_testout | wc -l"))
    elif shift<0:
        return int(p("(dd if=/dev/zero bs=%d count=1; cat /s/bpsk31_testin)>/s/bpsk31_testin0; cmp -l /s/bpsk31_testin0 /s/bpsk31_testout | wc -l"%-shift))
    elif shift>0:
        return int(p("(dd if=/dev/zero bs=%d count=1; cat /s/bpsk31_testout)>/s/bpsk31_testout0; cmp -l /s/bpsk31_testin /s/bpsk31_testout0 | wc -l"%shift))


lf=open("/s/output_results","w")

for snr in range(0,20,2):
    genfiles(snr)
    num_totalbits=getminsize()
    num_errors=None
    for shift in range(-5,5):
        curr_num_errors = mkdiff(shift)
        if not num_errors or (num_errors and num_errors > curr_num_errors):
            num_errors = curr_num_errors
    lf.write("%d; %d; %d; %d\n" %(snr, num_errors, num_totalbits, num_errors/float(num_totalbits)))
