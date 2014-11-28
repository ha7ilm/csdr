#ifdef USE_RPI_FFT

//I had an idea to use the GPU based FFT on Raspberry Pi boards.
//It would speed up filtering and spectrum display.
//However, this feature is not implemented yet. I've just started to work on it.

FFT_PLAN_T* make_fft_c(int size, complexf* input, complexf* output, int forward, int benchmark)
{
	int hmailbox = mbox_open();
	int returned = gpu_fft_prepare(hmailbox, log2N, (forward)?GPU_FFT_FWD:GPU_FFT_REV, jobs, &fft);
}

void fft_execute(FFT_PLAN_T* plan)
{
}

void fft_destroy(FFT_PLAN_T* plan)
{
}

#endif
