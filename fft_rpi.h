#pragma once

#ifdef USE_RPI_FFT
#include "mailbox.h"
#include "gpu_fft.h"

#define FFT_PLAN_T (struct GPU_FFT)

//reference: https://github.com/raspberrypi/userland/blob/master/host_applications/linux/apps/hello_pi/hello_fft/hello_fft.c

#endif
