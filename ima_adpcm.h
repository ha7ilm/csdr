#pragma once

#ifdef USE_IMA_ADPCM

typedef struct ImaState {
   int index;    // Index into step size table
   int previousValue; // Most recent sample value
} ima_adpcm_state_t;

ima_adpcm_state_t encode_ima_adpcm_i16_u8(short* input, unsigned char* output, int input_length, ima_adpcm_state_t state);
ima_adpcm_state_t decode_ima_adpcm_u8_i16(unsigned char* input, short* output, int input_length, ima_adpcm_state_t state);

#endif
