#!/bin/sh
# Run this script on a Raspberry Pi 2, while running test_shift_remote.grc on your PC. 
# It allows you to debug the NEON-accelerated version of specific DSP algorithms on the target hardware.
TEMPSCRIPT="/tmp/test_shift_remote_exec.sh"
echo '#!/bin/sh\ncsdr shift_addfast_cc -0.1' > $TEMPSCRIPT
cat $TEMPSCRIPT
chmod +x $TEMPSCRIPT
ncat -vvl 5321 -e $TEMPSCRIPT
rm $TEMPSCRIPT
